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

  # Rotate among four strong start-pattern templates to couple colors across rounds
  template_bank = [
    (2, 6, 10, 14),  # classic KT
    (1, 5, 9, 13),   # left-shifted
    (3, 7, 11, 15),  # right-shifted
    (4, 8, 12, 16),  # stretched-right
  ]

  T = [(0, 1)]
  for round_idx in range(6):
    starts = template_bank[round_idx % 4]
    lo = min(l for l, r in T)
    hi = max(r for l, r in T)
    delta = hi - lo if hi > lo else 1
    # build translated blocks
    blocks = []
    for start in starts:
      blocks.append([(delta * start + l - lo, delta * start + r - lo) for l, r in T])
    S = []
    # Interleave on even rounds, sequential on odd rounds
    if round_idx % 2 == 0:
      maxlen = max(len(b) for b in blocks)
      for i in range(maxlen):
        for blk in blocks:
          if i < len(blk):
            S.append(blk[i])
    else:
      for blk in blocks:
        S.extend(blk)
    # connectors based on the active starts
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

  # Micro-phase A: insert three long-range caps near the tail to boost FF without raising omega.
  lo = min(l for l, r in T)
  hi = max(r for l, r in T)
  delta = hi - lo if hi > lo else 1
  def cap_at(a_frac, b_frac):
    L = lo + max(1, int(round(a_frac * delta)))
    R = lo + max(1, int(round(b_frac * delta)))
    if R <= L:
      R = L + 1
    return (L, R)
  caps = [cap_at(0.08, 0.60), cap_at(0.25, 0.75), cap_at(0.75, 0.92)]
  if len(T) < CAP - 8:
    out = list(T)
    for i, cap in enumerate(caps):
      pos = len(out) - (i * 2 + 1)
      if pos < 0:
        out.append(cap)
      else:
        out.insert(pos, cap)
    T = out

  # Stage 2: Two delta2-driven micro extension rounds (thin sampling; four-start translations).
  # Skip if near capacity.
  if len(T) < CAP - 16:
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
      # Use half-scale to avoid expanding the global span too aggressively
      delta2 = max(1, G // 2)

      # Thin seed: bounded and deterministic size to respect budget
      per_block_target = max(8, min(64, budget // 12))
      U = thin_seed(current_T, per_block_target)
      if not U:
        return []

      # Build four translated blocks using the same 4-start template,
      # with deterministic parity-based internal reversal
      blocks = []
      ulo = min(l for l, r in U)
      s0, s1, s2, s3 = template_bank[0]
      for s in (s0, s1, s2, s3):
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
      connectors = [
        ((s0 - 1) * delta2 + glo, (s1 - 1) * delta2 + glo),  # left cap
        ((s2 + 2) * delta2 + glo, (s3 + 2) * delta2 + glo),  # right cap
        ((s0 + 2) * delta2 + glo, (s2 - 1) * delta2 + glo),  # cross 1
        ((s1 + 2) * delta2 + glo, (s3 - 1) * delta2 + glo),  # cross 2
        ((s0 + 3) * delta2 + glo, (s3 + 3) * delta2 + glo),  # cross3 (long range)
      ]
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

    # Micro-phase C: deterministic four-window augmentation (CAP-guarded)
    if len(T) < CAP - 8:
      def window_phase(current_T, budget, iter_id=0):
        if budget <= 0 or not current_T:
          return []
        glo = min(l for l, r in current_T)
        ghi = max(r for l, r in current_T)
        G = max(1, ghi - glo)

        # Thin, evenly spaced seed to limit omega growth
        seed_sz = max(8, min(32, budget // 4))
        U = thin_seed(current_T, seed_sz)
        if not U:
          return []
        ulo = min(l for l, r in U)

        # Four interior windows to avoid end-stack cliques
        window_fracs = [(0.15, 0.25), (0.32, 0.42), (0.58, 0.68), (0.82, 0.92)]
        blocks = []
        for (fa, fb) in window_fracs:
          win_lo = glo + int(round(fa * G))
          base = win_lo - ulo
          block = [(l + base, r + base) for (l, r) in U]
          # small deterministic reversal to break symmetry
          if (((int(fa * 100)) // 5 + iter_id) % 2) == 1:
            block = list(reversed(block))
          blocks.append(block)

        # Interleave the window blocks
        micro = []
        maxlen = max(len(b) for b in blocks)
        order = list(range(len(blocks)))
        if iter_id % 2 == 1:
          order.reverse()
        for i in range(maxlen):
          for idx in order:
            blk = blocks[idx]
            if i < len(blk):
              micro.append(blk[i])

        # Fractional-span connectors across windows
        micro_connectors = [
          (glo + int(round(0.08 * G)), glo + int(round(0.30 * G))),  # left cap
          (glo + int(round(0.60 * G)), glo + int(round(0.92 * G))),  # right cap
          (glo + int(round(0.26 * G)), glo + int(round(0.56 * G))),  # cross 1
          (glo + int(round(0.44 * G)), glo + int(round(0.78 * G))),  # cross 2
        ]
        for a, b in micro_connectors:
          if b > a:
            micro.append((a, b))

        if len(micro) > budget:
          micro = micro[:budget]
        return micro

      remaining = CAP - len(T)
      microC = window_phase(T, budget=min(196, remaining), iter_id=2)
      if microC:
        room = CAP - len(T)
        if len(microC) > room:
          microC = microC[:room]
        T.extend(microC)

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