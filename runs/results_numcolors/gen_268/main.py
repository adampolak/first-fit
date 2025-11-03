# EVOLVE-BLOCK-START

def construct_intervals(seed_count=1):
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

  # Use rotating KT start-pattern templates for strong coupling and diversity
  template_bank = [
    (2, 6, 10, 14),
    (1, 5, 9, 13),
    (3, 7, 11, 15),
    (4, 8, 12, 16),
  ]
  spine_starts = template_bank[0]

  # Seed with multiple disjoint unit intervals if requested (new capability)
  if seed_count <= 1:
    T = [(0, 1)]
  else:
    step = 3
    T = [(i * step, i * step + 1) for i in range(seed_count)]

  # Stage 1: six KT rounds, alternating interleaving to raise FF pressure
  for ridx in range(6):
    # select rotating template from bank
    starts = template_bank[ridx % len(template_bank)]
    # Predict next size: size -> 4*size + 4; abort if it would exceed CAP
    nxt_size = 4 * len(T) + 4
    if nxt_size > CAP:
      break
    lo = min(l for l, r in T)
    hi = max(r for l, r in T)
    delta = hi - lo
    if delta <= 0:
      delta = 1
    # Build four translated blocks
    blocks = []
    for start in starts:
      base = start * delta - lo
      blocks.append([(l + base, r + base) for (l, r) in T])
    # Interleave on even rounds, sequential on odd rounds
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
    # Classic connectors (Figure 4 style)
    s0, s1, s2, s3 = starts
    connectors = [
      ((s0 - 1) * delta, (s1 - 1) * delta),  # left cap
      ((s2 + 2) * delta, (s3 + 2) * delta),  # right cap
      ((s0 + 2) * delta, (s2 - 1) * delta),  # cross 1
      ((s1 + 2) * delta, (s3 - 1) * delta),  # cross 2
    ]
    S.extend(connectors)
    T = S

  # If we are already near the cap, return the strong baseline.
  if len(T) >= CAP - 16:
    return T

  # Stage 2: Two delta2-driven micro extension rounds (thin sampling; four-start translations).
  # Goals:
  #  - Raise FF pressure via cross-scale coupling and interleaving parity,
  #  - Keep omega in check by using thin seeds and sparse caps/connectors,
  #  - Respect strict capacity guard.

  def thin_seed(current_T, max_seed):
    """Take a thin, evenly spaced sample of current_T of size <= max_seed."""
    n = len(current_T)
    if n == 0 or max_seed <= 0:
      return []
    step = max(1, n // max_seed)
    U = current_T[::step][:max_seed]
    return U

  def micro_round(current_T, round_id, budget):
    if budget <= 0 or not current_T:
      return []

    glo = min(l for l, r in current_T)
    ghi = max(r for l, r in current_T)
    G = max(1, ghi - glo)

    # Use delta2 at half-scale to avoid expanding the global span too aggressively.
    delta2 = max(1, G // 2)

    # Thin seed: bounded and deterministic size to respect budget
    # Keep micro blocks small: target at most ~ (budget//6) per block (four blocks + ~8 extras).
    per_block_target = max(8, min(64, budget // 12))
    U = thin_seed(current_T, per_block_target)

    if not U:
      return []

    # Build four translated blocks using the same 4-start template,
    # with deterministic parity-based interleaving policy:
    # - Rounds with round_id % 2 == 0: forward interleave
    # - Rounds with round_id % 2 == 1: reverse interleave
    blocks = []
    ulo = min(l for l, r in U)
    for s in spine_starts:
      base = s * delta2 - ulo
      block = [(l + base, r + base) for (l, r) in U]
      # Add a tiny deterministic internal reversal to break symmetry every other block
      if ((s // 2) % 2) == (round_id % 2):
        block = list(reversed(block))
      blocks.append(block)

    micro = []
    if round_id % 2 == 0:
      # Forward interleave
      maxlen = max(len(b) for b in blocks)
      for i in range(maxlen):
        for blk in blocks:
          if i < len(blk):
            micro.append(blk[i])
    else:
      # Reverse interleave (swap block order)
      maxlen = max(len(b) for b in blocks)
      for i in range(maxlen):
        for blk in reversed(blocks):
          if i < len(blk):
            micro.append(blk[i])

    # Deterministic connectors at delta2 scale, including a cross3 extension
    s0, s1, s2, s3 = spine_starts
    connectors = [
      ((s0 - 1) * delta2 + glo, (s1 - 1) * delta2 + glo),  # left cap (localized)
      ((s2 + 2) * delta2 + glo, (s3 + 2) * delta2 + glo),  # right cap
      ((s0 + 2) * delta2 + glo, (s2 - 1) * delta2 + glo),  # cross 1
      ((s1 + 2) * delta2 + glo, (s3 - 1) * delta2 + glo),  # cross 2
      ((s0 + 3) * delta2 + glo, (s3 + 3) * delta2 + glo),  # cross3 (longer range)
    ]
    micro.extend(connectors)

    # Sparse micro caps: add three long but safe caps at half-scale
    cap1 = (glo + (delta2 // 4), glo + int(1.8 * delta2))
    cap2 = (glo + int(0.9 * delta2), glo + int(2.6 * delta2))
    mid = glo + G // 2
    cap3 = (mid - max(1, delta2 // 8), mid + max(1, delta2 // 8))
    for cap in (cap1, cap2, cap3):
      if cap[1] > cap[0]:
        micro.append(cap)

    # Enforce budget
    if len(micro) > budget:
      micro = micro[:budget]
    return micro

  # Execute two micro rounds with strong guards
  remaining = CAP - len(T)
  mr1 = micro_round(T, round_id=0, budget=max(0, remaining // 2))
  if mr1:
    room = CAP - len(T)
    if len(mr1) > room:
      mr1 = mr1[:room]
    T.extend(mr1)

  remaining = CAP - len(T)
  mr2 = micro_round(T, round_id=1, budget=max(0, remaining))
  if mr2:
    room = CAP - len(T)
    if len(mr2) > room:
      mr2 = mr2[:room]
    T.extend(mr2)

  return T

  # return [  # Figure 3, OPT=2, FF=4
  #   (2,3),
  #   (6,7),
  #   (10,11),
  #   (14,15),
  #   (1,5),
  #   (12,16),
  #   (4,9),
  #   (8,13),
  # ]

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()