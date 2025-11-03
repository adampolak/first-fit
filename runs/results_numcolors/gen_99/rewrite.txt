# EVOLVE-BLOCK-START

def construct_intervals(max_intervals=9800,
                        depth=6,
                        rotate_starts=True,
                        reverse_block_parity=False,
                        interleave_blocks=False,
                        micro_phase=False):
  """
  Rotated four-block Kierstead–Trotter style expansion with adaptive depth cap
  and an optional safeguarded micro-phase of sparse caps.

  Returns a list of (l, r) integer tuples (open intervals) in the order they
  should be presented to FirstFit.

  Parameters:
    max_intervals (int): hard cap on the total number of intervals produced.
    depth (int): nominal recursion depth (will be adaptively capped to respect max_intervals).
    rotate_starts (bool): rotate the four translated starts across rounds to disrupt repeating overlap patterns.
    reverse_block_parity (bool): reverse the order of T for every odd block within a round to mix color usage.
    interleave_blocks (bool): interleave blocks round-robin to maximize color mixing.
    micro_phase (bool): append a very small safeguarded micro-phase at the end with sparse long caps.

  Returns:
    intervals: list of tuples, each tuple (l, r) represents an open interval from l to r
  """

  # Capacity guard to keep total < max_intervals
  CAP = int(max_intervals)

  # Deterministic four-pattern cycle (four starts per round)
  start_patterns = [
    (2, 6, 10, 14),  # classic 4-wave
    (1, 5, 9, 13),   # left-shifted
    (3, 7, 11, 15),  # right-shifted
    (2, 4, 8, 12),   # compressed variant
  ]

  # Seed with a single unit interval; keep small to allow many rounds within budget
  T = [(0, 1)]

  # Predictive size accounting to cap the number of full rounds.
  # KT growth per round: size -> 4*size + 4 (four translated copies + four connectors)
  def round_next_size(sz):
    return 4 * sz + 4

  depth = max(0, int(depth))
  size = len(T)
  allowed = 0
  for _ in range(depth):
    nxt = round_next_size(size)
    if nxt > CAP:
      break
    size = nxt
    allowed += 1
  depth = allowed

  def apply_round(current_T, starts, do_interleave):
    # Compute span and delta (integers)
    lo = min(l for l, r in current_T)
    hi = max(r for l, r in current_T)
    delta = hi - lo
    if delta <= 0:
      delta = 1

    # Build translated blocks
    blocks = []
    for b_idx, s in enumerate(starts):
      # Optional reverse parity for odd-indexed block to disrupt symmetry
      src = current_T[::-1] if (reverse_block_parity and (b_idx % 2 == 1)) else current_T
      base = s * delta - lo
      block = [(l + base, r + base) for (l, r) in src]
      blocks.append(block)

    # Combine blocks: interleave or sequential
    S = []
    if do_interleave:
      maxlen = max(len(b) for b in blocks)
      order = list(range(len(blocks)))
      # small deterministic rotation by round parity for variety
      for i in range(maxlen):
        for idx in order:
          blk = blocks[idx]
          if i < len(blk):
            S.append(blk[i])
    else:
      for blk in blocks:
        S.extend(blk)

    # Classic KT connectors based on the first four starts
    s0, s1, s2, s3 = starts[:4]
    connectors = [
      ((s0 - 1) * delta, (s1 - 1) * delta),  # left cap
      ((s2 + 2) * delta, (s3 + 2) * delta),  # right cap
      ((s0 + 2) * delta, (s2 - 1) * delta),  # cross 1
      ((s1 + 2) * delta, (s3 - 1) * delta),  # cross 2
    ]
    S.extend(connectors)
    return S

  # Main KT rounds
  for round_idx in range(depth):
    if not T:
      break
    # Choose starts (rotate if requested)
    starts = start_patterns[round_idx % len(start_patterns)] if rotate_starts else start_patterns[0]
    # Apply round; by default keep sequential blocks which has shown strong performance
    T = apply_round(T, starts, do_interleave=bool(interleave_blocks))

    # Capacity guard mid-process
    if len(T) > CAP:
      T = T[:CAP]
      break

  # Optional safeguarded micro-phase: add a tiny number of long caps if requested
  if micro_phase and T and len(T) + 2 <= CAP:
    lo = min(l for l, r in T)
    hi = max(r for l, r in T)
    span = max(1, hi - lo)
    d2 = max(1, span // 3)
    # Two sparse caps away from densest core, inspired by the safer scaffold
    micro_caps = [
      (lo + 1 * d2, lo + 5 * d2),
      (hi - 6 * d2, hi - 2 * d2),
    ]
    # Insert near the end to influence late FF assignments without densifying cliques
    for i, itv in enumerate(micro_caps):
      pos = max(0, len(T) - 1 - i)
      T.insert(pos, itv)

  # Normalize to non-negative integers
  if not T:
    return []
  min_l = min(l for l, r in T)
  if min_l < 0:
    T = [(l - min_l, r - min_l) for (l, r) in T]

  # Ensure integer open intervals and monotone right > left
  intervals = []
  for (l, r) in T:
    li = int(l)
    ri = int(r)
    if ri <= li:
      ri = li + 1
    intervals.append((li, ri))

  # Final capacity trim
  if len(intervals) > CAP:
    intervals = intervals[:CAP]

  return intervals

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()