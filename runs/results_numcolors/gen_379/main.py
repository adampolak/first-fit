# EVOLVE-BLOCK-START

def construct_intervals(enable_alt_microphase=True):
  """
  Construct a sequence of intervals of real line,
  in the order in which they are presented to FirstFit,
  so that it maximizes the number of colors used by FirstFit
  divided by the maximum number of intervals that cover a single point.

  Deterministic KT-style backbone with two CAP-aware micro-phases and
  carefully chosen connectors, keeping omega small (target <= 10).

  Args:
    enable_alt_microphase (bool): enable a second, alternate micro-phase.

  Returns:
    intervals: list of tuples, each tuple (l, r) represents an open interval from l to r.
  """

  # Capacity guard to keep total intervals < 10000
  CAP = 9800

  # Rotating six-start templates for the KT spine (expanded bank for more mixing)
  template_bank = [
    (2, 6, 10, 14),  # classic KT
    (1, 5, 9, 13),   # left-shifted
    (3, 7, 11, 15),  # right-shifted
    (4, 8, 12, 16),  # stretched-right
    (5, 9, 13, 17),  # extra right-shift
    (6, 10, 14, 18), # extra classic-shift
  ]

  # Seed with a single unit interval
  T = [(0, 1)]

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

    # Interleaving policy: order derived from 'starts' to diversify per-round mixing
    S = []
    if do_interleave:
      # define order by sorting block indices by start values (or reverse)
      order = sorted(range(4), key=lambda i: starts[i], reverse=reverse_order)
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

    # Classic four connectors (Figure-4 style)
    s0, s1, s2, s3 = starts
    S.append(((s0 - 1) * delta, (s1 - 1) * delta))  # left cap
    S.append(((s2 + 2) * delta, (s3 + 2) * delta))  # right cap
    S.append(((s0 + 2) * delta, (s2 - 1) * delta))  # cross 1
    S.append(((s1 + 2) * delta, (s3 - 1) * delta))  # cross 2
    return S

  # Stage 1: KT spine with rotation and parity policies, capacity-guarded
  for ridx in range(8):
    # Predict next size: sz -> 4*sz + 4
    if 4 * len(T) + 4 > CAP:
      break
    starts = template_bank[ridx % len(template_bank)]
    # Interleave on even rounds, reverse block order on odd rounds
    do_inter = (ridx % 2 == 0)
    rev = (ridx % 2 == 1)
    T = _apply_round(T, starts, do_interleave=do_inter, reverse_order=rev)

  # Early exit if nearly at capacity (before micro-phases)
  if len(T) >= CAP - 8:
    return T

  # Spine densification: duplicate a thin seed of intervals with slight shrink to boost FF pressure
  room = CAP - len(T)
  if room > 16:
    # Sample a small seed of intervals uniformly from the KT spine
    stride = max(1, len(T) // 100)
    seed = [T[i] for i in range(0, len(T), stride)][:10]
    lo_d, hi_d, span_d = _span(T)
    eps = max(1, span_d // 1000)
    densify = []
    for (l, r) in seed:
      if r - l > 2 * eps:
        densify.append((l + eps, r - eps))
    # Insert densified intervals if capacity allows
    for iv in densify:
      if len(T) >= CAP:
        break
      T.append(iv)

  # Micro-phase A: add three long-range caps to boost FF pressure (capacity-safe)
  lo, hi, _ = _span(T)
  span = max(1, hi - lo)
  def cap_at(a_frac, b_frac):
    L = lo + max(1, int(round(a_frac * span)))
    R = lo + max(1, int(round(b_frac * span)))
    if R <= L:
      R = L + 1
    return (L, R)
  caps = [cap_at(0.08, 0.60), cap_at(0.25, 0.75), cap_at(0.75, 0.92)]
  for c in caps:
    if len(T) >= CAP:
      break
    T.append(c)

  # Micro-phase builder: thin-window replication with connectors
  def build_micro_delta_round(current_T, budget, alt=False):
    if not current_T or budget <= 8:
      return []

    glo = min(l for l, r in current_T)
    ghi = max(r for l, r in current_T)
    G = max(1, ghi - glo)

    # Thin, evenly spaced seed
    seed_sz = max(8, min(32, len(current_T) // 300))
    if alt:
      seed_sz = max(8, min(40, len(current_T) // 250))
    stride = max(1, len(current_T) // max(1, seed_sz))
    U = [current_T[i] for i in range(0, len(current_T), stride)][:seed_sz]
    if not U:
      return []

    ulo = min(l for l, r in U)

    # Four windows across the span (stay inside to control omega)
    if not alt:
      window_fracs = [(0.12, 0.22), (0.35, 0.45), (0.58, 0.68), (0.80, 0.90)]
    else:
      window_fracs = [(0.05, 0.15), (0.28, 0.38), (0.60, 0.70), (0.82, 0.92)]

    # Build translated micro-blocks aligned to these windows
    blocks = []
    for (fa, fb) in window_fracs:
      win_lo = glo + int(round(fa * G))
      base = win_lo - ulo
      block = [(l + base, r + base) for (l, r) in U]
      # Internal reversal on alternating pattern for alt-phase to diversify
      if alt and int(round(fa * 100)) // 5 % 2 == 1:
        block = list(reversed(block))
      blocks.append(block)

    # Interleave micro-blocks to maximize FF mixing
    micro = []
    maxlen = max(len(b) for b in blocks)
    if not alt:
      # forward interleave
      for i in range(maxlen):
        for blk in blocks:
          if i < len(blk):
            micro.append(blk[i])
    else:
      # reverse interleave
      for i in range(maxlen):
        for blk in reversed(blocks):
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

    # A long cross to tie distant colors (alt-phase only)
    if alt:
      a = glo + int(round(0.18 * G))
      b = glo + int(round(0.84 * G))
      if b > a:
        micro.append((a, b))

    # Trim to available budget
    if len(micro) > budget:
      micro = micro[:budget]
    return micro

  # Micro-phase B: CAP-aware first micro-round; reserve budget for later phases
  room = CAP - len(T)
  if room > 8:
    budgetB = max(16, int(room * 0.55))
    micro = build_micro_delta_round(T, budgetB, alt=False)
    if micro:
      if len(micro) > budgetB:
        micro = micro[:budgetB]
      T.extend(micro)

  # Micro-phase C: alternate window family + long cross (guarded, with explicit budget)
  room = CAP - len(T)
  if enable_alt_microphase and room > 8:
    budgetC = max(24, min(room - 8, int(room * 0.35), 512))
    micro_alt = build_micro_delta_round(T, budgetC, alt=True)
    if micro_alt:
      if len(micro_alt) > budgetC:
        micro_alt = micro_alt[:budgetC]
      T.extend(micro_alt)

  # Micro-phase D: sparse ladder of disjoint long-range connectors (off-center windows)
  # to couple many active colors without stacking at the spine's peak clique region.
  room = CAP - len(T)
  if room > 6:
    glo = min(l for l, r in T)
    ghi = max(r for l, r in T)
    G = max(1, ghi - glo)
    # Six disjoint windows away from the central peak to avoid omega spikes
    win_fracs = [(0.10, 0.18), (0.28, 0.36), (0.46, 0.52),
                 (0.62, 0.68), (0.76, 0.82), (0.88, 0.94)]
    ladder = []
    for (fa, fb) in win_fracs:
      a = glo + int(round(fa * G))
      b = glo + int(round(fb * G))
      if b > a:
        ladder.append((a, b))
    if ladder:
      ladder = ladder[:min(len(ladder), room)]
      T.extend(ladder)

  # Micro-phase E: final fractional connectors; add two short edge-biased segments
  room = CAP - len(T)
  if room > 0:
    glo_, hi_, delta_ = _span(T)
    connectors_q = [
      (glo_ + delta_ // 4, hi_ - delta_ // 4),
      (glo_ + delta_ // 2 - delta_ // 8, glo_ + delta_ // 2 + delta_ // 8),
      (glo_ + 3 * delta_ // 4 - 1, hi_ - 1),
      (glo_ + delta_ // 12, glo_ + delta_ // 5),
      (glo_ + 9 * delta_ // 10, glo_ + 11 * delta_ // 12),
    ]
    for c in connectors_q:
      if len(T) < CAP:
        T.append(c)
  # Final capacity trim
  if len(T) > CAP: T = T[:CAP]

  return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()