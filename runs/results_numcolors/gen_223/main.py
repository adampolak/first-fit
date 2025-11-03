# EVOLVE-BLOCK-START

def construct_intervals(seed_count=1):
  """
  Construct a sequence of intervals on the real line, in the order presented to FirstFit,
  to maximize FF colors divided by the clique number.
  Returns:
    intervals: list of (l, r) tuples (open intervals).
  """

  # Capacity guard to keep total intervals < 10000
  CAP = 9800

  # Rotating four-start templates for the KT spine
  template_bank = [
    (2, 6, 10, 14),  # classic KT
    (1, 5, 9, 13),   # left-shifted
    (3, 7, 11, 15),  # right-shifted
    (4, 8, 12, 16),  # stretched-right
  ]

  # Seed with a single unit interval; multi-seed tends to inflate omega too early
  T = [(0, 1)] if seed_count <= 1 else [(0, 1)]

  # Helpers
  def _span(current_T):
    lo = min(l for l, r in current_T)
    hi = max(r for l, r in current_T)
    delta = hi - lo
    if delta <= 0:
      delta = 1
    return lo, hi, delta

  def _apply_round(current_T, starts, do_interleave=False, reverse_order=False):
    lo, hi, delta = _span(current_T)

    # Build four translated blocks
    blocks = []
    for s in starts:
      base = s * delta - lo
      block = [(l + base, r + base) for (l, r) in current_T]
      blocks.append(block)

    # Interleaving policy (round-robin across blocks)
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
    return S

  # Stage 1: KT spine with rotation and parity policies, capacity-guarded
  for ridx in range(6):
    # Predict next size: sz -> 4*sz + 4
    if 4 * len(T) + 4 > CAP:
      break
    starts = template_bank[ridx % len(template_bank)]
    # Interleave on even rounds, reverse block order on odd rounds
    do_inter = (ridx % 2 == 0)
    rev = (ridx % 2 == 1)
    T = _apply_round(T, starts, do_interleave=do_inter, reverse_order=rev)

  # Early exit if nearly at capacity
  if len(T) >= CAP - 8:
    return T

  # Micro-phase A: add three long-range caps to boost FF pressure (capacity-safe)
  lo, hi, _ = _span(T)
  span = max(1, hi - lo)
  def cap_at(a_frac, b_frac):
    L = lo + max(1, int(round(a_frac * span)))
    R = lo + max(1, int(round(b_frac * span)))
    if R <= L:
      R = L + 1
    return (L, R)
  # Add a fourth long cross cap to couple distant colors while keeping omega tame
  caps = [cap_at(0.08, 0.60), cap_at(0.25, 0.75), cap_at(0.75, 0.92), cap_at(0.18, 0.84)]
  for c in caps:
    if len(T) >= CAP:
      break
    T.append(c)

  # Micro-phase B: thin delta2-style fractional-window round(s), capacity-guarded
  def build_micro_delta_round(current_T, budget):
    if not current_T or budget <= 8:
      return []

    glo = min(l for l, r in current_T)
    ghi = max(r for l, r in current_T)
    G = max(1, ghi - glo)

    # Thin, evenly spaced seed
    seed_sz = max(8, min(32, len(current_T) // 300))
    stride = max(1, len(current_T) // max(1, seed_sz))
    U = [current_T[i] for i in range(0, len(current_T), stride)][:seed_sz]
    if not U:
      return []

    ulo = min(l for l, r in U)

    # Four windows across the span (stay inside to control omega)
    window_fracs = [(0.12, 0.22), (0.35, 0.45), (0.58, 0.68), (0.80, 0.90)]

    # Build translated micro-blocks aligned to these windows
    blocks = []
    for (fa, fb) in window_fracs:
      win_lo = glo + int(round(fa * G))
      base = win_lo - ulo
      block = [(l + base, r + base) for (l, r) in U]
      blocks.append(block)

    # Interleave micro-blocks to maximize FF mixing
    micro = []
    maxlen = max(len(b) for b in blocks)
    for i in range(maxlen):
      for blk in blocks:
        if i < len(blk):
          micro.append(blk[i])

    # Fractional-span connectors at the micro scale
    micro_connectors = [
      (glo + int(round(0.08 * G)), glo + int(round(0.30 * G))),  # left cap
      (glo + int(round(0.60 * G)), glo + int(round(0.92 * G))),  # right cap
      (glo + int(round(0.26 * G)), glo + int(round(0.56 * G))),  # cross 1
      (glo + int(round(0.44 * G)), glo + int(round(0.78 * G))),  # cross 2
    ]
    for a, b in micro_connectors:
      if b > a:
        micro.append((a, b))

    # Trim to available budget
    if len(micro) > budget:
      micro = micro[:budget]
    return micro

  # One safeguarded micro-round (empirically strong baseline)
  room = CAP - len(T)
  if room > 8:
    micro = build_micro_delta_round(T, room)
    if micro:
      if len(micro) > room:
        micro = micro[:room]
      T.extend(micro)

  # Second guarded micro-phase with alternate window family (fixed, well inside span)
  def build_micro_delta_round_alt(current_T, budget):
    if not current_T or budget <= 8:
      return []

    glo = min(l for l, r in current_T)
    ghi = max(r for l, r in current_T)
    G = max(1, ghi - glo)

    # Slightly larger thin seed to strengthen coupling while bounded
    seed_sz = max(8, min(40, len(current_T) // 250))
    stride = max(1, len(current_T) // max(1, seed_sz))
    U = [current_T[i] for i in range(0, len(current_T), stride)][:seed_sz]
    if not U:
      return []

    ulo = min(l for l, r in U)

    # Alternate window set; chosen to avoid boundaries and control omega
    window_fracs = [(0.05, 0.15), (0.28, 0.38), (0.60, 0.70), (0.82, 0.92)]

    # Build translated micro-blocks aligned to these windows with internal parity reversals
    blocks = []
    for (fa, fb) in window_fracs:
      win_lo = glo + int(round(fa * G))
      base = win_lo - ulo
      block = [(l + base, r + base) for (l, r) in U]
      # Break symmetry on alternating windows deterministically
      if (int(round(fa * 100)) // 5) % 2 == 1:
        block = list(reversed(block))
      blocks.append(block)

    # Reverse-order interleaving to diversify cross-block mixing
    micro = []
    maxlen = max(len(b) for b in blocks)
    for i in range(maxlen):
      for blk in reversed(blocks):
        if i < len(blk):
          micro.append(blk[i])

    # Fractional-span connectors including a long cross to tie distant colors
    micro_connectors = [
      (glo + int(round(0.08 * G)), glo + int(round(0.30 * G))),
      (glo + int(round(0.60 * G)), glo + int(round(0.92 * G))),
      (glo + int(round(0.26 * G)), glo + int(round(0.56 * G))),
      (glo + int(round(0.44 * G)), glo + int(round(0.78 * G))),
      (glo + int(round(0.18 * G)), glo + int(round(0.84 * G))),
    ]
    for a, b in micro_connectors:
      if b > a:
        micro.append((a, b))

    if len(micro) > budget:
      micro = micro[:budget]
    return micro

  # Execute alternate micro-phase if capacity allows
  room = CAP - len(T)
  if room > 8:
    micro_alt = build_micro_delta_round_alt(T, room)
    if micro_alt:
      if len(micro_alt) > room:
        micro_alt = micro_alt[:room]
      T.extend(micro_alt)

  # Final capacity trim
  if len(T) > CAP:
    T = T[:CAP]

  return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()