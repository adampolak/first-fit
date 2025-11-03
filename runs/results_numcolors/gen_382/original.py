# EVOLVE-BLOCK-START

def _span_delta(T):
  lo = min(l for l, _ in T)
  hi = max(r for _, r in T)
  delta = hi - lo
  if delta <= 0:
    delta = 1
  return lo, hi, delta

def _append_connectors(S, starts, delta, add_cross4=False):
  s0, s1, s2, s3 = starts
  # Classic four connectors (safe and strong)
  S.append(((s0 - 1) * delta, (s1 - 1) * delta))  # left cap
  S.append(((s2 + 2) * delta, (s3 + 2) * delta))  # right cap
  S.append(((s0 + 2) * delta, (s2 - 1) * delta))  # cross 1
  S.append(((s1 + 2) * delta, (s3 - 1) * delta))  # cross 2
  # Optional long-range cross (conservative; use only near the end)
  if add_cross4:
    S.append(((s0 + 4) * delta, (s3 + 4) * delta))

def _apply_round(current_T, starts, do_interleave=False, reverse_order=False, add_cross4=False):
  lo, hi, delta = _span_delta(current_T)

  # Build translated blocks from current_T
  base_offsets = [s * delta - lo for s in starts]
  blocks = [[(l + base, r + base) for (l, r) in current_T] for base in base_offsets]

  # Assemble S with optional interleaving and reverse order to mix colors
  S = []
  if do_interleave:
    maxlen = max(len(b) for b in blocks)
    order = list(range(4))
    if reverse_order:
      order.reverse()
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

  # Classic connectors (+ optional long-range cross)
  _append_connectors(S, starts, delta, add_cross4=add_cross4)
  return S

def construct_intervals(seed_count=1, phase2_iters=2, enable_parabolic=False):
  """
  Construct a sequence of intervals of real line,
  in the order in which they are presented to FirstFit,
  so that it maximizes the number of colors used by FirstFit
  divided by the maximum number of intervals that cover a single point.

  Hybrid: six-round KT spine with rotated starts and parity-based assembly,
  followed by a near-tail cap injection and thin delta-driven micro rounds.
  """
  CAP = 9800

  # Rotated start templates (empirically strong)
  template_bank = [
    (2, 6, 10, 14),  # classic KT
    (1, 5, 9, 13),   # left-shifted
    (3, 7, 11, 15),  # right-shifted
    (4, 8, 12, 16),  # stretched-right
  ]

  # Seed with a single unit interval to avoid early omega inflation
  T = [(0, 1)]

  # Stage 1: up to six KT rounds, alternating interleaving with reversed block order on odd rounds
  # Add a conservative long-range cross only on the final feasible round to enhance coupling.
  max_rounds = 6
  for ridx in range(max_rounds):
    starts = template_bank[ridx % len(template_bank)]
    nxt_size = 4 * len(T) + 4
    if nxt_size > CAP:
      break
    do_interleave = (ridx % 2 == 0)
    reverse_order = (ridx % 2 == 1)
    # Enable long cross4 on the final round only if we still leave modest headroom
    # This mirrors the safer policy from the inspiration program.
    will_have_room = (4 * nxt_size + 4) <= CAP  # rough lookahead
    add_cross4 = (ridx == max_rounds - 1) and (not do_interleave) and will_have_room
    T = _apply_round(T, starts, do_interleave=do_interleave, reverse_order=reverse_order, add_cross4=add_cross4)
    if len(T) >= CAP:
      T = T[:CAP]
      return T

  # If near capacity, return strong baseline
  if len(T) >= CAP - 16:
    return T

  # Micro-phase A: inject sparse long caps near the tail to couple many active colors
  def _insert_near_tail(seq, intervals):
    out = list(seq)
    for i, iv in enumerate(intervals):
      pos = len(out) - (i * 2 + 1)
      if pos < 0:
        out.append(iv)
      else:
        out.insert(pos, iv)
    return out

  lo, hi, delta = _span_delta(T)
  def cap_at(a_frac, b_frac):
    span = max(1, hi - lo)
    L = lo + max(1, int(round(a_frac * span)))
    R = lo + max(1, int(round(b_frac * span)))
    if R <= L:
      R = L + 1
    return (L, R)

  tail_caps = [
    cap_at(0.08, 0.60),
    cap_at(0.25, 0.75),
    cap_at(0.75, 0.92),
  ]
  room = CAP - len(T)
  if room > 0:
    T = _insert_near_tail(T, tail_caps[:room])

  if len(T) >= CAP - 16:
    return T

  # Micro-phase B: delta-driven micro rounds using a thin evenly spaced seed
  def build_micro_delta_round(current_T, budget, iter_id=0, alt=False):
    if not current_T or budget <= 8:
      return []

    glo = min(l for l, r in current_T)
    ghi = max(r for l, r in current_T)
    G = max(1, ghi - glo)

    # Thin seed (slightly larger for stronger mixing but bounded)
    seed_sz = max(8, min(40, len(current_T) // 250))
    stride = max(1, len(current_T) // max(1, seed_sz))
    U = [current_T[i] for i in range(0, len(current_T), stride)][:seed_sz]
    if not U:
      return []

    ulo = min(l for l, r in U)

    # Window families (primary with slight shifts; alternate is distinct and wider)
    if not alt:
      shift = (iter_id % 3) * 0.02
      window_fracs = [
        (0.12 + shift, 0.22 + shift),
        (0.35 + shift, 0.45 + shift),
        (0.58 + shift, 0.68 + shift),
        (0.80 + shift, 0.90 + shift),
      ]
      window_fracs = [
        (max(0.05, min(0.90, a)), max(0.10, min(0.95, b)))
        for (a, b) in window_fracs
      ]
    else:
      window_fracs = [(0.05, 0.15), (0.28, 0.38), (0.60, 0.70), (0.82, 0.92)]

    # Build translated micro-blocks aligned to windows
    blocks = []
    for (fa, fb) in window_fracs:
      win_lo = glo + int(round(fa * G))
      base = win_lo - ulo
      block = [(l + base, r + base) for (l, r) in U]
      tag = iter_id if not alt else (iter_id + 1)
      if ((int(round(fa * 100)) // 5) + tag) % 2 == 1:
        block = list(reversed(block))
      blocks.append(block)

    # Interleave micro-blocks with tag-dependent order
    micro = []
    maxlen = max(len(b) for b in blocks)
    order = list(range(len(blocks)))
    tag = iter_id if not alt else (iter_id + 1)
    if tag % 2 == 1:
      order.reverse()
    for i in range(maxlen):
      for idx in order:
        blk = blocks[idx]
        if i < len(blk):
          micro.append(blk[i])

    # Connectors at fractional span scale (alternate adds a long cross)
    micro_connectors = [
      (glo + int(round(0.08 * G)), glo + int(round(0.30 * G))),  # left cap
      (glo + int(round(0.60 * G)), glo + int(round(0.92 * G))),  # right cap
      (glo + int(round(0.26 * G)), glo + int(round(0.56 * G))),  # cross 1
      (glo + int(round(0.44 * G)), glo + int(round(0.78 * G))),  # cross 2
    ]
    if alt:
      micro_connectors.append((glo + int(round(0.18 * G)), glo + int(round(0.84 * G))))
    for a, b in micro_connectors:
      if b > a:
        micro.append((a, b))

    if len(micro) > budget:
      micro = micro[:budget]
    return micro

  # Execute micro rounds; phase2_iters serves as the requested iterations (capped at 2)
  steps = min(max(0, int(phase2_iters)), 2)
  for iter_id in range(steps):
    room = CAP - len(T)
    if room <= 8:
      break
    micro = build_micro_delta_round(T, room, iter_id=iter_id, alt=False)
    if not micro:
      break
    if len(micro) > room:
      micro = micro[:room]
    T.extend(micro)

  # Alternate micro family to capture missed interactions
  room = CAP - len(T)
  if room > 8:
    microB = build_micro_delta_round(T, room, iter_id=steps, alt=True)
    if microB:
      if len(microB) > room:
        microB = microB[:room]
      T.extend(microB)

  # Optional lightweight parabolic caps
  if enable_parabolic:
    room = CAP - len(T)
    if room > 0:
      lo, hi, delta = _span_delta(T)
      caps = [
        (lo + max(1, int(round(0.15 * delta))), lo + max(1, int(round(0.60 * delta)))),
        (lo + max(1, int(round(0.25 * delta))), lo + max(1, int(round(0.80 * delta)))),
        (lo + max(1, int(round(0.55 * delta))), lo + max(1, int(round(0.95 * delta)))),
      ]
      add = [c for c in caps if c[1] > c[0]][:room]
      T.extend(add)

  return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()