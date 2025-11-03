# EVOLVE-BLOCK-START

def construct_intervals(max_intervals=9800,
                        depth=6,
                        rotate_starts=True,
                        reverse_block_parity=False,
                        interleave_blocks=True,
                        micro_phase=True):
  """
  Efficient, deterministic KT-style interval spine with optional micro-phase.
  Returns a list of (l, r) tuples representing open intervals in the order
  intended for input to FirstFit.

  Maintains compatibility with the original signature and output format.
  """

  # Cap on total intervals to respect the provided budget
  CAP = int(max_intervals)

  # Four-block start patterns (classic KT spine plus three alternates)
  start_patterns = [
    (2, 6, 10, 14),  # A
    (1, 5, 9, 13),   # B
    (3, 7, 11, 15),  # C
    (2, 4, 8, 12),   # D
  ]

  # Seed with a single small interval to maximize round count within CAP
  T = [(0, 1)]

  # Helper: round growth estimate to cap depth adaptively
  def round_next_size(sz):
    # Each round adds 4 copies plus 4 connectors: 4*sz + 4
    return 4 * sz + 4

  # Compute how many rounds can be safely performed within CAP
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

  # Build one KT-round given the current spine starts
  def apply_round(current_T, starts, do_interleave):
    # Span and delta
    lo = min(l for l, r in current_T)
    hi = max(r for l, r in current_T)
    delta = hi - lo
    if delta <= 0:
      delta = 1

    # Translate current_T by each start; optionally reverse parity per block
    blocks = []
    for idx, s in enumerate(starts):
      base = s * delta - lo
      if reverse_block_parity and (idx % 2 == 1):
        # reverse source T for this block to break symmetry
        src = current_T[::-1]
      else:
        src = current_T
      block = [(l + base, r + base) for (l, r) in src]
      blocks.append(block)

    # Merge blocks: interleave if requested, otherwise concatenate
    S = []
    if do_interleave:
      maxlen = max((len(b) for b in blocks), default=0)
      for i in range(maxlen):
        for blk in blocks:
          if i < len(blk):
            S.append(blk[i])
    else:
      for blk in blocks:
        S.extend(blk)

    # Classic connectors derived from the first four starts
    s0, s1, s2, s3 = starts[:4]
    connectors = [
      ((s0 - 1) * delta, (s1 - 1) * delta),  # left cap
      ((s2 + 2) * delta, (s3 + 2) * delta),  # right cap
      ((s0 + 2) * delta, (s2 - 1) * delta),  # cross 1
      ((s1 + 2) * delta, (s3 - 1) * delta),  # cross 2
      ((s0 + 3) * delta, (s3 + 3) * delta),  # cross 3
    ]
    S.extend(connectors)
    return S

  # Main KT rounds with optional rotation and interleaving
  for round_idx in range(depth):
    if not T:
      break
    starts = start_patterns[round_idx % len(start_patterns)] if rotate_starts else start_patterns[0]
    T = apply_round(T, starts, do_interleave=bool(interleave_blocks))
    if len(T) > CAP:
      T = T[:CAP]
      break

  # Optional micro-phase: small, sparse long caps added near the end
  if micro_phase and T and len(T) + 3 <= CAP:
    lo = min(l for l, r in T)
    hi = max(r for l, r in T)
    span = max(1, hi - lo)
    d2 = max(1, span // 3)
    center = lo + span // 2
    width = max(1, span // 16)

    micro_caps = [
      (lo + 1 * d2, lo + 5 * d2),
      (hi - 6 * d2, hi - 2 * d2),
      (center - width, center + width),  # center cap
    ]
    # Append near end to influence late FF without dense cliques
    T.extend(micro_caps)

  # Normalize to non-negative integers, ensure valid intervals
  if not T:
      return []

  min_l = min(l for l, r in T)
  if min_l < 0:
    T = [(l - min_l, r - min_l) for (l, r) in T]

  intervals = []
  for (l, r) in T:
    li = int(l)
    ri = int(r)
    if ri <= li:
      ri = li + 1
    intervals.append((li, ri))

  if len(intervals) > CAP:
    intervals = intervals[:CAP]

  return intervals

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()