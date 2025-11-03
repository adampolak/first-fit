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

  # Use the classic fixed start-pattern to preserve strong KT coupling
  starts = (2, 6, 10, 14)
  T = [(0, 1)]
  for round_idx in range(6):
    lo = min(l for l, r in T)
    hi = max(r for l, r in T)
    delta = hi - lo
    # build translated blocks in sequence (no interleaving)
    blocks = []
    for start in starts:
      blocks.append([(delta * start + l - lo, delta * start + r - lo) for l, r in T])
    S = []
    # keep sequential order every round to stabilize omega and reproduce KT backbone
    for blk in blocks:
      S.extend(blk)
    # connectors based on the fixed starts
    s0, s1, s2, s3 = starts
    connectors = [
      (delta * (s0 - 1), delta * (s1 - 1)),  # left cap
      (delta * (s2 + 2), delta * (s3 + 2)),  # right cap
      (delta * (s0 + 2), delta * (s2 - 1)),  # cross 1
      (delta * (s1 + 2), delta * (s3 - 1)),  # cross 2
    ]
    for a, b in connectors:
      S.append((a, b))
    T = S
  # Stage 2: Two delta2-driven micro extension rounds (thin sampling; four-start translations).
  # Skip if near capacity.
  if len(T) < CAP - 16:
    def thin_seed(current_T, max_seed):
      n = len(current_T)
      if n == 0 or max_seed <= 0:
        return []
      step = max(1, n // max_seed)
      return current_T[::step][:max_seed]

    def micro_round(current_T, round_id, budget, starts_param, scale_den=2, include_cross4=False):
      if budget <= 0 or not current_T:
        return []
      glo = min(l for l, r in current_T)
      ghi = max(r for l, r in current_T)
      G = max(1, ghi - glo)
      # Use reduced-scale to avoid expanding the global span too aggressively
      delta2 = max(1, G // max(1, int(scale_den)))

      # Thin seed: bounded and deterministic size to respect budget
      per_block_target = max(8, min(64, budget // 12))
      U = thin_seed(current_T, per_block_target)
      if not U:
        return []

      # Build translated blocks using the provided start template,
      # with deterministic parity-based internal reversal
      blocks = []
      ulo = min(l for l, r in U)
      starts_local = starts_param
      for s in starts_local:
        base = s * delta2 - ulo
        block = [(l + base, r + base) for (l, r) in U]
        if ((s // 2) % 2) == (round_id % 2):
          block = list(reversed(block))
        blocks.append(block)

      micro = []
      # Parity-based interleaving across blocks
      maxlen = max(len(b) for b in blocks)
      if round_id % 2 == 0:
        for i in range(maxlen):
          for blk in blocks:
            if i < len(blk):
              micro.append(blk[i])
      else:
        for i in range(maxlen):
          for blk in reversed(blocks):
            if i < len(blk):
              micro.append(blk[i])

      # Deterministic connectors at delta2 scale, including a longer cross3
      s0, s1, s2, s3 = starts_local[0], starts_local[1], starts_local[2], starts_local[3]
      connectors = [
        ((s0 - 1) * delta2 + glo, (s1 - 1) * delta2 + glo),  # left cap
        ((s2 + 2) * delta2 + glo, (s3 + 2) * delta2 + glo),  # right cap
        ((s0 + 2) * delta2 + glo, (s2 - 1) * delta2 + glo),  # cross 1
        ((s1 + 2) * delta2 + glo, (s3 - 1) * delta2 + glo),  # cross 2
        ((s0 + 3) * delta2 + glo, (s3 + 3) * delta2 + glo),  # cross3 (long range)
      ]
      if include_cross4:
        connectors.append(((s0 + 1) * delta2 + glo, (s3 - 1) * delta2 + glo))  # cross4 (very long)
      micro.extend(connectors)

      # Sparse micro caps to press FirstFit locally without raising omega much
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

    remaining = CAP - len(T)
    mr1 = micro_round(T, round_id=0, budget=max(0, remaining // 3), starts_param=starts, scale_den=2, include_cross4=False)
    if mr1:
      room = CAP - len(T)
      if len(mr1) > room:
        mr1 = mr1[:room]
      T.extend(mr1)

    remaining = CAP - len(T)
    mr2 = micro_round(T, round_id=1, budget=max(0, remaining // 2), starts_param=starts, scale_den=2, include_cross4=False)
    if mr2:
      room = CAP - len(T)
      if len(mr2) > room:
        mr2 = mr2[:room]
      T.extend(mr2)

    # Second micro-phase with distinct starts and smaller scale to amplify cross-scale coupling
    remaining = CAP - len(T)
    alt_starts = (1, 5, 9, 13)
    mr3 = micro_round(T, round_id=2, budget=max(0, remaining), starts_param=alt_starts, scale_den=3, include_cross4=True)
    if mr3:
      room = CAP - len(T)
      if len(mr3) > room:
        mr3 = mr3[:room]
      T.extend(mr3)

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