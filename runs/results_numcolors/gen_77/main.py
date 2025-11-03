# EVOLVE-BLOCK-START

def construct_intervals(rounds=3):
  """
  Construct a sequence of open intervals presented to FirstFit.
  Strategy:
    - Seed with a low-omega spine of disjoint unit intervals.
    - For each round:
        * Translate the current set by multiple deterministic start patterns.
        * Interleave block order with parity variation to break color reuse.
        * Add deterministic connectors ("caps") to couple colors across blocks.
        * Run a bounded micro second-phase using delta2 to expand cross-scale ties.
    - Normalize to integers and cap to a safe maximum size.

  Args:
    rounds (int): number of macro expansion rounds.

  Returns:
    List[Tuple[int, int]]: open intervals (l, r) in presentation order.
  """
  # Hard cap to ensure tractability and evaluator bounds.
  MAX_INTERVALS = 9800

  # Deterministic bank of start-pattern templates (eight-round cycle).
  # Each entry is a list of start offsets used (scaled by delta).
  start_bank = [
    [2, 6, 10, 14],        # A
    [1, 5, 9, 13],         # B
    [3, 7, 11, 15],        # C
    [2, 4, 8, 12],         # D (compressed)
    [6, 10, 14, 18],       # E (right-extended)
    [5, 9, 13, 17],        # F (left-extended)
    [7, 11, 15, 19],       # G (farther right)
    [4, 8, 12, 16],        # H (mid-dense)
  ]

  # Seed: four disjoint unit intervals as a low-omega spine
  T = [(0, 1), (2, 3), (4, 5), (6, 7)]

  # Utility: ensure budget before extending S
  def can_add(cur_size, add_count):
    return (cur_size + add_count) <= MAX_INTERVALS

  # Main rounds
  for round_idx in range(max(0, int(rounds))):
    if not T:
      break

    lo = min(l for l, r in T)
    hi = max(r for l, r in T)
    delta = hi - lo
    if delta <= 0:
      delta = 1

    # Choose a start pattern based on round and also tentatively extend density on later rounds
    starts = list(start_bank[round_idx % len(start_bank)])
    # When budget allows and round >= 1, extend with two extra offsets to enrich interactions
    if round_idx >= 1:
      ext = [starts[-1] + 4, starts[-1] + 8]
      starts.extend(ext)

    # Budget-aware trimming of starts to avoid blow-up
    if len(T) * len(starts) > MAX_INTERVALS // 2:
      # Trim the pattern deterministically
      keep = max(3, min(len(starts), MAX_INTERVALS // max(1, (2 * len(T)))))
      starts = starts[:keep]

    # Build translated blocks with parity reversal for odd blocks
    blocks = []
    for b_idx, s in enumerate(starts):
      src = T if ((b_idx % 2) == 0) else list(reversed(T))
      base = s * delta - lo
      # Fast path when budget gets tight
      blk = [(l + base, r + base) for (l, r) in src]
      blocks.append(blk)

    # Interleave blocks round-robin with deterministic ordering
    S = []
    if blocks:
      maxlen = max(len(b) for b in blocks)
      order = list(range(len(blocks)))
      # Deterministic rotation and direction flip every two rounds
      if (round_idx // 2) % 2 == 1:
        order = order[::-1]
      rot = (2 * round_idx) % len(order)
      order = order[rot:] + order[:rot]

      for i in range(maxlen):
        # Alternate inner traversal direction by round
        idxs = order if (round_idx % 2 == 0) else list(reversed(order))
        for idx in idxs:
          blk = blocks[idx]
          if i < len(blk):
            if can_add(len(S), 1):
              S.append(blk[i])
            else:
              break

    # Deterministic connectors ("caps") to couple colors across blocks without large cliques
    if starts:
      s0 = starts[0]
      s1 = starts[1] if len(starts) > 1 else s0 + 4
      s2 = starts[-2] if len(starts) > 2 else s0 + 8
      s3 = starts[-1] if len(starts) > 1 else s0 + 12

      connectors = [
        ((s0 - 1) * delta, (s1 - 1) * delta),  # left cap
        ((s2 + 2) * delta, (s3 + 2) * delta),  # right cap
        ((s0 + 2) * delta, (s2 - 1) * delta),  # cross 1
        ((s1 + 2) * delta, (s3 - 1) * delta),  # cross 2
      ]
      # Add a middle-range long connector to engage more colors at low clique
      if len(starts) >= 4:
        connectors.append(((s0 + 1.5) * delta, (s3 - 1.5) * delta))

      for a, b in connectors:
        if a < b and can_add(len(S), 1):
          S.append((a, b))

    # Secondary delta2-driven micro-phase
    # Use small sub-blocks to create cross-scale coupling without exploding size
    if S:
      loS = min(l for l, r in S)
      hiS = max(r for l, r in S)
      span = hiS - loS
      delta2 = int(max(1, span // 6))  # coarser than delta to reduce overlap inflation

      # Choose a bounded slice of S for micro-expansion
      slice_len = min(64, max(16, len(S) // 64))
      S_head = S[:slice_len]

      # Two micro rounds with alternating small start patterns
      micro_patterns = ([1, 3], [2, 5])
      for micro_idx, mp in enumerate(micro_patterns):
        additions = []
        for ms in mp:
          base2 = ms * delta2 - loS
          # Add a sparse sub-block translated by base2
          for (l, r) in S_head:
            if can_add(len(S) + len(additions), 1):
              additions.append((l + base2, r + base2))
            else:
              break
        # A pair of micro-connectors per micro round
        mc_a = (loS + 1 * delta2, loS + 3 * delta2)
        mc_b = (hiS - 4 * delta2, hiS - 2 * delta2)
        for mc in (mc_a, mc_b):
          if mc[0] < mc[1] and can_add(len(S) + len(additions), 1):
            additions.append(mc)
        # Append micro additions
        S.extend(additions)

    # If nearing budget, stop expanding further rounds
    if len(S) >= MAX_INTERVALS:
      T = S[:MAX_INTERVALS]
      break

    T = S

  # Normalize to non-negative integers and ensure open intervals (ri > li)
  if not T:
    return []

  min_l = min(l for l, r in T)
  if min_l < 0:
    T = [(l - min_l, r - min_l) for l, r in T]

  intervals = []
  intervals_reserve = min(len(T), MAX_INTERVALS)
  for i in range(intervals_reserve):
    l, r = T[i]
    li = int(round(l))
    ri = int(round(r))
    if ri <= li:
      ri = li + 1
    intervals.append((li, ri))

  return intervals

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()