# EVOLVE-BLOCK-START

def construct_intervals(max_intervals=9800,
                        target_depth=6,
                        rotate_starts=True,
                        reverse_block_parity=True,
                        interleave_blocks=True,
                        seed_size=3,
                        micro_phase=True):
  """
  Crossover implementation combining multi-block Kierstead-style expansion
  with seeds, rotated starts, interleaving, block-parity reversal, per-round
  lightweight blockers, and an optional micro-phase of long caps.

  Returns a list of (l, r) integer tuples (open intervals) in FF presentation order.

  Parameters:
    max_intervals (int): hard cap on total intervals returned.
    target_depth (int): nominal recursion depth.
    rotate_starts (bool): rotate start patterns across rounds.
    reverse_block_parity (bool): reverse every other block's order in a round.
    interleave_blocks (bool): interleave copies round-robin each round.
    seed_size (int): number of small disjoint seed intervals to start with.
    micro_phase (bool): enable final sprinkling of long caps.

  The construction is deterministic and aims to pressure FirstFit while keeping
  the clique number moderate.
  """
  # Bank of start-pattern templates (we will rotate among them if requested)
  start_patterns = [
    [2, 6, 10, 14],     # classic
    [1, 5, 9, 13],      # left-shift
    [3, 7, 11, 15],     # right-shift
    [2, 4, 8, 12],      # compressed
  ]

  # Seed: several small disjoint unit intervals to diversify initial colors
  # keep separation 2 so they don't form cliques
  T = [(i * 3, i * 3 + 1) for i in range(max(1, int(seed_size)))]

  # Estimate per-round growth conservatively: copies_per_round = len(starts),
  # connectors = 4, blockers = min(4, copies). We'll cap depth such that
  # size after depth rounds <= max_intervals.
  copies_per_round = len(start_patterns[0])  # normally 4
  connectors_per_round = 4
  blockers_per_round = copies_per_round  # we may add one small blocker per copy

  # Compute adaptive depth
  depth = max(0, int(target_depth))
  size = len(T)
  allowed = 0
  while allowed < depth:
    next_size = size * copies_per_round + connectors_per_round + blockers_per_round
    if next_size > max_intervals:
      break
    size = next_size
    allowed += 1
  depth = allowed

  # Helper to append small blockers (short intervals inside block spans) without
  # increasing global clique significantly. They are placed localized to each block.
  def make_blockers(delta, lo, starts):
    blockers = []
    # Use short intervals of length 1 positioned inside each translated block span.
    # Offset them so they don't overlap across different blocks (reducing omega growth).
    for i, s in enumerate(starts):
      base = s * delta - lo
      # offset within block: place at quarter-point of the base shift
      offset = (i % 3) + 1  # small integer shift varying per block
      l = base + offset
      r = l + 1
      blockers.append((l, r))
    return blockers

  # Main recursive/iterative expansion rounds
  for round_idx in range(depth):
    lo = min(l for l, r in T)
    hi = max(r for l, r in T)
    delta = hi - lo
    if delta <= 0:
      delta = 1

    # Choose start tuple for this round (rotate to break patterns if requested)
    starts = start_patterns[round_idx % len(start_patterns)] if rotate_starts else start_patterns[0]

    # Create translated copies (with optional reversal per block to mix color order)
    blocks = []
    for b_idx, s in enumerate(starts):
      block_src = T[::-1] if (reverse_block_parity and (b_idx % 2 == 1)) else T
      base = s * delta - lo
      block = [(l + base, r + base) for (l, r) in block_src]
      blocks.append(block)

    # Interleave blocks round-robin (or append block-by-block)
    S = []
    if interleave_blocks:
      maxlen = max(len(b) for b in blocks)
      # rotate order each round for variety
      order = list(range(len(blocks)))
      if round_idx % 2 == 1:
        order = order[::-1]
      krot = round_idx % len(order)
      order = order[krot:] + order[:krot]
      for i in range(maxlen):
        for idx in order:
          blk = blocks[idx]
          if i < len(blk):
            S.append(blk[i])
    else:
      for blk in blocks:
        S.extend(blk)

    # Append lightweight blockers to consume early small colors in each block
    blockers = make_blockers(delta, lo, starts)
    # Insert blockers interleaved to have them appear before connectors
    for b in blockers:
      S.append(b)

    # Deterministic connectors derived from the first four starts to propagate color chains
    s0, s1, s2, s3 = starts[:4]
    connectors = [
      ((s0 - 1) * delta, (s1 - 1) * delta),  # left cap
      ((s2 + 2) * delta, (s3 + 2) * delta),  # right cap
      ((s0 + 2) * delta, (s2 - 1) * delta),  # cross 1
      ((s1 + 2) * delta, (s3 - 1) * delta),  # cross 2
    ]
    for a, b in connectors:
      S.append((a, b))

    T = S

    # Safety cap: if T exceeds max_intervals, truncate early
    if len(T) > max_intervals:
      T = T[:max_intervals]
      break

  # Optional micro-phase: sprinkle a few sparse long-range caps near the end
  if micro_phase:
    if T:
      lo = min(l for l, r in T)
      hi = max(r for l, r in T)
      delta = max(1, hi - lo)
      d2 = max(1, delta // 5)
      micro = [
        (lo + 1 * d2, lo + 5 * d2),
        (hi - 6 * d2, hi - 2 * d2),
        (lo + 3 * d2, lo + 8 * d2),
      ]
      # Insert micro gadgets near the end of the sequence to touch many colors
      insert_base = max(0, len(T) - 3)
      for i, it in enumerate(micro):
        pos = insert_base + i
        if pos >= len(T):
          T.append(it)
        else:
          T.insert(pos, it)

  # Normalize to non-negative integers and ensure integer open intervals of positive length
  if not T:
    return []
  min_l = min(l for l, r in T)
  if min_l < 0:
    T = [(l - min_l, r - min_l) for l, r in T]

  intervals = []
  for (l, r) in T:
    # Round to nearest integer; ensure length >=1
    li = int(round(l))
    ri = int(round(r))
    if ri <= li:
      ri = li + 1
    intervals.append((li, ri))

  # Final size cap
  if len(intervals) > max_intervals:
    intervals = intervals[:max_intervals]

  return intervals

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()