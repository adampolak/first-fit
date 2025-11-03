# EVOLVE-BLOCK-START

def construct_intervals(rounds=6,
                        rotate_starts=True,
                        reverse_block_parity=True,
                        interleave_blocks=True,
                        phase2_iters=1,
                        cross4_enabled=True):
  """
  Rotating-template KT spine with dual fractional-window micro phases and guarded cross4.
  Returns intervals (l, r) in FirstFit presentation order.
  """

  # Capacity guard to keep total intervals < 10000
  CAP = 9800

  # Start template(s)
  spine_starts = (2, 6, 10, 14)
  template_bank = [
    (2, 6, 10, 14),  # classic KT
    (1, 5, 9, 13),   # left-shifted
    (3, 7, 11, 15),  # right-shifted
    (4, 8, 12, 16),  # stretched-right
  ]

  # Seed with one unit interval
  T = [(0, 1)]

  # Helpers
  def _span(current_T):
    lo = min(l for l, r in current_T)
    hi = max(r for l, r in current_T)
    delta = hi - lo
    if delta <= 0:
      delta = 1
    return lo, hi, delta

  def _apply_round(current_T, starts, do_interleave=False, reverse_order=False, add_cross4=False):
    lo, hi, delta = _span(current_T)

    # Build four translated blocks
    blocks = []
    for s in starts:
      base = s * delta - lo
      block = [(l + base, r + base) for (l, r) in current_T]
      blocks.append(block)

    # Interleaving policy
    S = []
    if do_interleave:
      order = list(range(4))
      if reverse_order:
        order.reverse()
      maxlen = max(len(b) for b in blocks)
      for i in range(maxlen):
        for idx in order:
          blk = blocks[idx]
          if i < len(blk):
            S.append(blk[i])
    else:
      if reverse_order:
        blocks = list(reversed(blocks))
      for blk in blocks:
        S.extend(blk)

    # Classic four connectors
    s0, s1, s2, s3 = starts
    S.append(((s0 - 1) * delta, (s1 - 1) * delta))  # left cap
    S.append(((s2 + 2) * delta, (s3 + 2) * delta))  # right cap
    S.append(((s0 + 2) * delta, (s2 - 1) * delta))  # cross 1
    S.append(((s1 + 2) * delta, (s3 - 1) * delta))  # cross 2

    # Optional long-range cross4 connector (single-scale, conservative)
    if add_cross4:
      S.append(((s0 + 4) * delta, (s3 + 4) * delta))
    return S

  # Stage 1: KT spine with rotation and parity policies, capacity-guarded.
  # We enable a single cross4 connector on the last executed round if cross4_enabled is True.
  executed_rounds = 0
  pending_rounds = max(0, int(rounds))
  for ridx in range(pending_rounds):
    # Predict next size: sz -> 4*sz + 4
    if 4 * len(T) + 4 > CAP:
      break
    starts = template_bank[ridx % len(template_bank)] if rotate_starts else spine_starts
    do_inter = bool(interleave_blocks and (ridx % 2 == 0))
    rev = bool(reverse_block_parity and (ridx % 2 == 1))
    # Defer decision about cross4; add only on the last round actually executed.
    add_c4 = False
    # Apply round
    T = _apply_round(T, starts, do_interleave=do_inter, reverse_order=rev, add_cross4=add_c4)
    executed_rounds += 1

  # Retrofit: add one cross4 on the last round scale by appending it now (safer than per-round accumulation)
  if cross4_enabled and executed_rounds > 0:
    lo, hi, delta = _span(T)
    # Approximate the last used starts; we align with the deterministic base template used last.
    last_starts = template_bank[(executed_rounds - 1) % len(template_bank)] if rotate_starts else spine_starts
    s0, _, _, s3 = last_starts
    # Add a single long connector at the outermost scale (conservative)
    cross4_iv = ((s0 + 4) * delta, (s3 + 4) * delta)
    T.append(cross4_iv)

  # Early exit if nearly at capacity
  if len(T) >= CAP - 64:
    return T[:CAP]

  # Micro-phase A: add three long-range caps to boost FF pressure (capacity-safe)
  lo, hi, _ = _span(T)
  G = max(1, hi - lo)

  def cap_at(a_frac, b_frac):
    L = lo + max(1, int(round(a_frac * G)))
    R = lo + max(1, int(round(b_frac * G)))
    if R <= L:
      R = L + 1
    return (L, R)

  caps = [cap_at(0.08, 0.60), cap_at(0.25, 0.75), cap_at(0.75, 0.92)]
  for c in caps:
    if len(T) >= CAP:
      break
    T.append(c)

  # Micro-phase builder with parameterized windows and deterministic sampling offset
  def build_micro_delta_round(current_T, budget, windows, interleave_reverse=False, seed_offset=0, add_cross4_local=False):
    if not current_T or budget <= 8:
      return []

    glo = min(l for l, r in current_T)
    ghi = max(r for l, r in current_T)
    Gspan = max(1, ghi - glo)

    # Thin, evenly spaced seed with a deterministic offset
    seed_sz = max(8, min(32, len(current_T) // 300))
    stride = max(1, len(current_T) // max(1, seed_sz))
    # Offset the start index deterministically to diversify micro phases
    start_idx = (seed_offset % stride)
    U = [current_T[i] for i in range(start_idx, len(current_T), stride)][:seed_sz]
    if not U:
      return []

    ulo = min(l for l, r in U)

    # Build translated micro-blocks aligned to these windows
    blocks = []
    for (fa, fb) in windows:
      # Clamp windows inside (0.04, 0.96) to avoid accidental boundary cliques
      fa = max(0.04, min(0.92, fa))
      fb = max(fa + 0.04, min(0.96, fb))
      win_lo = glo + int(round(fa * Gspan))
      base = win_lo - ulo
      block = [(l + base, r + base) for (l, r) in U]
      blocks.append(block)

    # Interleave micro-blocks to maximize FF mixing with parity
    micro = []
    maxlen = max(len(b) for b in blocks)
    order = list(range(len(blocks)))
    if interleave_reverse:
      order.reverse()
    for i in range(maxlen):
      for idx in order:
        blk = blocks[idx]
        if i < len(blk):
          micro.append(blk[i])

    # Fractional-span connectors at the micro scale
    micro_connectors = [
      (glo + int(round(0.08 * Gspan)), glo + int(round(0.30 * Gspan))),  # left cap
      (glo + int(round(0.60 * Gspan)), glo + int(round(0.92 * Gspan))),  # right cap
      (glo + int(round(0.26 * Gspan)), glo + int(round(0.56 * Gspan))),  # cross 1
      (glo + int(round(0.44 * Gspan)), glo + int(round(0.78 * Gspan))),  # cross 2
    ]
    for a, b in micro_connectors:
      if b > a:
        micro.append((a, b))

    # Optional long-range cross4-like connector at micro scale (single instance)
    if add_cross4_local:
      a = glo + int(round(0.12 * Gspan))
      b = glo + int(round(0.94 * Gspan))
      if b > a:
        micro.append((a, b))

    # Trim to available budget
    if len(micro) > budget:
      micro = micro[:budget]
    return micro

  # Micro-phase B1: primary delta2 windows
  windows1 = [(0.12, 0.22), (0.35, 0.45), (0.58, 0.68), (0.80, 0.90)]
  room = CAP - len(T)
  if room > 8:
    micro1 = build_micro_delta_round(
      T, room, windows1,
      interleave_reverse=False,
      seed_offset=1,
      add_cross4_local=False
    )
    if micro1:
      if len(micro1) > room:
        micro1 = micro1[:room]
      T.extend(micro1)

  # Respect requested phase2_iters but guard to at most two micro-rounds.
  steps = min(max(0, int(phase2_iters)), 2)
  # If the user asked for only one micro-iteration, we already ran it; otherwise add the second variant.
  if steps >= 2 and len(T) < CAP - 8:
    windows2 = [(0.05, 0.15), (0.28, 0.38), (0.60, 0.70), (0.82, 0.92)]
    room = CAP - len(T)
    if room > 8:
      micro2 = build_micro_delta_round(
        T, room, windows2,
        interleave_reverse=True,
        seed_offset=3,
        add_cross4_local=True if cross4_enabled else False
      )
      if micro2:
        if len(micro2) > room:
          micro2 = micro2[:room]
        T.extend(micro2)

  # Final capacity trim
  if len(T) > CAP:
    T = T[:CAP]

  return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()