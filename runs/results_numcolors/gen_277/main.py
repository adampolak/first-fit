# EVOLVE-BLOCK-START

def construct_intervals(enable_alt_microphase=False):
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

  # Helper: conservative estimate of next round size (4*sz + connectors)
  def est_next_size(sz, add_cross4):
    return 4 * sz + (6 if add_cross4 else 4)

  for round_idx in range(6):
    # CAP-aware decision: try to enable cross4, but disable if it would overflow CAP
    add_cross4 = True
    if est_next_size(len(T), add_cross4) > CAP:
      add_cross4 = False
    if est_next_size(len(T), add_cross4) > CAP:
      break

    lo = min(l for l, r in T)
    hi = max(r for l, r in T)
    delta = hi - lo if hi > lo else 1

    # Spine density multiplier: alternate to increase cross-scale coupling on even rounds
    K = 2 if (round_idx % 2 == 0) else 1

    # build translated blocks
    blocks = []
    for start in starts:
      base = (start * K) * delta - lo
      blocks.append([(l + base, r + base) for l, r in T])

    S = []
    # parity interleaving: even rounds interleave blocks (round-robin), odd rounds sequential
    if round_idx % 2 == 0:
      maxlen = max(len(b) for b in blocks) if blocks else 0
      for i in range(maxlen):
        for blk in blocks:
          if i < len(blk):
            S.append(blk[i])
    else:
      for blk in blocks:
        S.extend(blk)

    # connectors based on the fixed starts and scaled by K
    s0, s1, s2, s3 = starts
    connectors = [
      ((s0 - 1) * K * delta, (s1 - 1) * K * delta),  # left cap
      ((s2 + 2) * K * delta, (s3 + 2) * K * delta),  # right cap
      ((s0 + 2) * K * delta, (s2 - 1) * K * delta),  # cross 1
      ((s1 + 2) * K * delta, (s3 - 1) * K * delta),  # cross 2
    ]
    for a, b in connectors:
      S.append((a, b))

    # Optional long-range cross4 connectors on safe rounds to transmit color pressure
    if add_cross4 and (round_idx % 2 == 0):
      cross4 = [
        ((s0 + 4) * K * delta, (s3 + 4) * K * delta),
        ((s1 + 4) * K * delta, (s2 + 4) * K * delta),
      ]
      for a, b in cross4:
        if b > a:
          S.append((a, b))

    # Small deterministic tie-pins to occupy low colors early and break symmetry
    eps = max(1, int(delta // 256))
    tie_positions = [lo + max(1, delta // 8), lo + max(1, delta // 2), lo + max(1, (3 * delta) // 4)]
    for pos in tie_positions:
      S.append((pos, pos + eps))

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

    def micro_round(current_T, round_id, budget, windows=None, long_cross=False):
      if budget <= 0 or not current_T:
        return []
      glo = min(l for l, r in current_T)
      ghi = max(r for l, r in current_T)
      G = max(1, ghi - glo)
      # choose delta2 to vary scale between micro rounds
      delta2 = max(1, G // (3 if (round_id % 2 == 0) else 2))

      # Thin seed: bounded and deterministic size to respect budget
      per_block_target = max(8, min(64, budget // 12))
      U = thin_seed(current_T, per_block_target)
      if not U:
        return []

      # Build four translated blocks using the same 4-start template,
      # with deterministic parity-based internal reversal
      blocks = []
      ulo = min(l for l, r in U)
      for s in starts:
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

      # Deterministic connectors at delta2 scale
      s0, s1, s2, s3 = starts
      connectors = [
        ((s0 - 1) * delta2 + glo, (s1 - 1) * delta2 + glo),  # left cap
        ((s2 + 2) * delta2 + glo, (s3 + 2) * delta2 + glo),  # right cap
        ((s0 + 2) * delta2 + glo, (s2 - 1) * delta2 + glo),  # cross 1
        ((s1 + 2) * delta2 + glo, (s3 - 1) * delta2 + glo),  # cross 2
      ]
      if long_cross:
        connectors.append((glo + int(round(0.18 * G)), glo + int(round(0.84 * G))))
      micro.extend(connectors)

      # Sparse micro caps to press FirstFit locally without raising omega much
      cap1 = (glo + (delta2 // 4), glo + int(1.8 * delta2))
      cap2 = (glo + int(0.9 * delta2), glo + int(2.6 * delta2))
      mid = glo + G // 2
      cap3 = (mid - max(1, delta2 // 8), mid + max(1, delta2 // 8))
      for cap in (cap1, cap2, cap3):
        if cap[1] > cap[0]:
          micro.append(cap)

      # Optional window-based micro blocks for added pressure
      if windows:
        for (fa, fb) in windows:
          win_lo = glo + int(round(fa * G))
          base = win_lo - ulo
          block = [(l + base, r + base) for (l, r) in U]
          micro.extend(block)

      # Enforce budget
      if len(micro) > budget:
        micro = micro[:budget]
      return micro

    remaining = CAP - len(T)
    # Micro-phase A: default windows with a long cross connector (if budget allows)
    window_fracs_A = [(0.12, 0.22), (0.35, 0.45), (0.58, 0.68), (0.80, 0.90)]
    mr1 = micro_round(T, round_id=0, budget=max(0, remaining // 2), windows=window_fracs_A, long_cross=True)
    if mr1:
      room = CAP - len(T)
      if len(mr1) > room:
        mr1 = mr1[:room]
      T.extend(mr1)

    remaining = CAP - len(T)
    # Micro-phase B: alternate windows (enabled by flag) to diversify coupling
    window_fracs_B = [(0.06, 0.14), (0.28, 0.38), (0.54, 0.64), (0.78, 0.88)]
    chosen_windows = window_fracs_B if enable_alt_microphase else window_fracs_A
    mr2 = micro_round(T, round_id=1, budget=max(0, remaining), windows=chosen_windows, long_cross=False)
    if mr2:
      room = CAP - len(T)
      if len(mr2) > room:
        mr2 = mr2[:room]
      T.extend(mr2)

  # Final integer conversion (keep values modest and strictly increasing endpoints)
  out = []
  for (l, r) in T[:CAP]:
    li = int(round(l))
    ri = int(round(r))
    if ri <= li:
      ri = li + 1
    out.append((li, ri))
  return out

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