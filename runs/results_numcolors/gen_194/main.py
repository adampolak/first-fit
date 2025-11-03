# EVOLVE-BLOCK-START

def construct_intervals():
  """
  Construct a sequence of intervals of real line,
  in the order in which they are presented to FirstFit,
  so that it maximizes the number of colors used by FirstFit
  divided by the maximum number of intervals that cover a single point

  The initial implementation uses the construction from
  Figure 4 in https://arxiv.org/abs/1506.00192

  Returns:
    intervals: list of tuples, each tuple (l, r) represents an open interval from l to r
  """

  # Capacity guard to keep total intervals < 10000
  CAP = 9800

  # Stage 1: Deterministic KT-style spine with rotating start patterns to diversify block interactions
  templates = [
      (2, 6, 10, 14),
      (1, 5, 9, 13),
      (3, 7, 11, 15),
      (4, 8, 12, 16),
  ]

  # Seed: keep a single unit interval; multi-seed tends to inflate omega too early.
  T = [(0, 1)]

  # Helper to insert caps near the tail of the presentation order (press late FF choices).
  def _insert_near_tail(seq, intervals):
    out = list(seq)
    for i, iv in enumerate(intervals):
      pos = len(out) - (i * 2 + 1)
      if pos < 0:
        out.append(iv)
      else:
        out.insert(pos, iv)
    return out

  # Perform up to six spine rounds unless capacity would be exceeded (it won't for this growth).
  for ridx in range(6):
    # Compute span and delta
    lo = min(l for l, r in T)
    hi = max(r for l, r in T)
    delta = hi - lo
    if delta <= 0:
      delta = 1

    # Rotate start pattern and build four translated blocks
    starts = templates[ridx % len(templates)]
    blocks = []
    for s in starts:
      base = s * delta - lo
      block = [(l + base, r + base) for (l, r) in T]
      blocks.append(block)

    # Interleave blocks on even rounds; keep sequential on odd rounds
    S = []
    if ridx % 2 == 0:
      maxlen = max(len(b) for b in blocks)
      for i in range(maxlen):
        for blk in blocks:
          if i < len(blk):
            S.append(blk[i])
    else:
      for blk in blocks:
        S.extend(blk)

    # Classic connectors (Figure 4 style), now using the rotated pattern
    s0, s1, s2, s3 = starts
    connectors = [
      ((s0 - 1) * delta, (s1 - 1) * delta),  # left cap
      ((s2 + 2) * delta, (s3 + 2) * delta),  # right cap
      ((s0 + 2) * delta, (s2 - 1) * delta),  # cross 1
      ((s1 + 2) * delta, (s3 - 1) * delta),  # cross 2
    ]

    # Add a cautious long-range cross4 connector on the final spine round only.
    # This couples colors across widely separated translated blocks without creating dense local cliques.
    if ridx == 5:
      connectors.append(((s0 + 4) * delta, (s3 + 4) * delta))

    S.extend(connectors)
    T = S

    if len(T) > CAP:
      break

  # Stage 2: delta2-driven micro extension rounds (thin sampling; four-start translations).
  # Goals:
  #  - Raise FF pressure via cross-scale coupling and interleaving parity,
  #  - Keep omega in check by using thin seeds and sparse caps/connectors,
  #  - Respect strict capacity guard and keep OPT <= 10 by avoiding dense local overlaps.
  remaining = CAP - len(T)
  if remaining <= 16:
    return T

  def thin_seed(current_T, max_seed):
    n = len(current_T)
    if n == 0 or max_seed <= 0:
      return []
    step = max(1, n // max_seed)
    return current_T[::step][:max_seed]

  def micro_round(current_T, round_id, budget):
    if budget <= 0 or not current_T:
      return []

    glo = min(l for l, r in current_T)
    ghi = max(r for l, r in current_T)
    G = max(1, ghi - glo)

    # Use a slightly smaller micro-scale (denser packing inside the span).
    delta2 = max(1, G // 3)

    # Thin seed: bounded and deterministic size to respect budget.
    per_block_target = max(8, min(96, max(12, budget // 10)))
    U = thin_seed(current_T, per_block_target)
    if not U:
      return []

    # Alternate window families across micro rounds for richer interactions
    if round_id % 3 == 0:
      window_fracs = [(0.12, 0.22), (0.35, 0.45), (0.58, 0.68), (0.80, 0.90)]
    elif round_id % 3 == 1:
      window_fracs = [(0.05, 0.15), (0.28, 0.38), (0.60, 0.70), (0.82, 0.92)]
    else:
      # denser central packing round
      window_fracs = [(0.10, 0.20), (0.30, 0.40), (0.50, 0.65), (0.70, 0.85)]

    # Build four translated blocks using a fixed 4-start template to stabilize omega.
    starts = (2, 6, 10, 14)
    blocks = []
    ulo = min(l for l, r in U)
    for i, s in enumerate(starts):
      # Place block aligned to a chosen window to diversify offsets
      win = window_fracs[i % len(window_fracs)]
      win_lo_frac = win[0]
      win_lo = glo + int(round(win_lo_frac * G))
      base = win_lo - ulo + s * delta2 - (s * delta2 // 4)
      block = [(l + base, r + base) for (l, r) in U]
      # Add a tiny deterministic internal reversal to break symmetry every other block
      if ((s // 2) % 2) == (round_id % 2):
        block = list(reversed(block))
      blocks.append(block)

    micro = []
    # Interleave with parity-dependent order to mix colors
    if round_id % 2 == 0:
      maxlen = max(len(b) for b in blocks)
      for i in range(maxlen):
        for blk in blocks:
          if i < len(blk):
            micro.append(blk[i])
    else:
      maxlen = max(len(b) for b in blocks)
      for i in range(maxlen):
        for blk in reversed(blocks):
          if i < len(blk):
            micro.append(blk[i])

    # Deterministic connectors at delta2 scale. Add longer-range cross connectors on later micro rounds.
    s0, s1, s2, s3 = starts
    connectors = [
      ((s0 - 1) * delta2 + glo, (s1 - 1) * delta2 + glo),  # left cap (localized)
      ((s2 + 2) * delta2 + glo, (s3 + 2) * delta2 + glo),  # right cap
      ((s0 + 2) * delta2 + glo, (s2 - 1) * delta2 + glo),  # cross 1
      ((s1 + 2) * delta2 + glo, (s3 - 1) * delta2 + glo),  # cross 2
    ]
    # Add extra long cross connectors only in later micro rounds to avoid early omega spikes.
    if round_id >= 1:
      connectors.append(((s0 + 3) * delta2 + glo, (s3 + 3) * delta2 + glo))  # cross3
    if round_id >= 2:
      connectors.append(((s1 + 4) * delta2 + glo, (s2 + 4) * delta2 + glo))  # cross4

    micro.extend(connectors)

    # Sparse micro caps: add three long but safe caps at micro scale to press many active colors.
    cap1 = (glo + max(1, delta2 // 6), glo + int(1.6 * delta2))
    cap2 = (glo + int(0.8 * delta2), glo + int(2.2 * delta2))
    mid = glo + G // 2
    cap3 = (mid - max(1, delta2 // 10), mid + max(1, delta2 // 10))
    for cap in (cap1, cap2, cap3):
      if cap[1] > cap[0]:
        micro.append(cap)

    # Enforce budget
    if len(micro) > budget:
      micro = micro[:budget]
    return micro

  # Execute up to three micro rounds (capacity-guarded)
  steps = 3
  for iter_id in range(steps):
    room = CAP - len(T)
    if room <= 8:
      break
    # allocate budget prudently: gradually increase per-round allowance
    per_budget = max(0, room // (steps - iter_id))
    micro = micro_round(T, round_id=iter_id, budget=per_budget)
    if not micro:
      continue
    if len(micro) > room:
      micro = micro[:room]
    T.extend(micro)

  # After micro rounds, insert a few near-tail long caps to press late FirstFit assignments
  lo = min(l for l, r in T)
  hi = max(r for l, r in T)
  delta = max(1, hi - lo)
  tail_caps = [
    (lo + max(1, int(0.07 * delta)), lo + max(1, int(0.58 * delta))),
    (lo + max(1, int(0.22 * delta)), lo + max(1, int(0.78 * delta))),
  ]
  room = CAP - len(T)
  if room > 0:
    T = _insert_near_tail(T, tail_caps[:room])

  # Ensure capacity
  if len(T) > CAP:
    T = T[:CAP]

  return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()