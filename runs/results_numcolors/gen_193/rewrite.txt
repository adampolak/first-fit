# EVOLVE-BLOCK-START

def construct_intervals(rounds=6,
                        rotate_starts=True,
                        reverse_block_parity=True,
                        interleave_blocks=True,
                        phase2_iters=1):
  """
  Deterministic KT-style spine with an improved, capacity-aware micro-phase.

  Parameters (same signature as original):
    rounds (int): main expansion depth; adaptive to CAP.
    rotate_starts (bool): rotate among templates when True.
    reverse_block_parity (bool): flip block source order for odd-indexed blocks.
    interleave_blocks (bool): whether to interleave blocks during spine rounds.
    phase2_iters (int): number of micro-phase iterations (safeguarded).

  Returns:
    intervals: list of (l, r) integer tuples, open intervals, in FF presentation order.
  """

  # Hard capacity guard to keep the total count < 10000
  CAP = 9800

  # Template bank / start patterns for rotation
  template_bank = [
    (2, 6, 10, 14),  # classic KT
    (1, 5, 9, 13),   # left biased
    (3, 7, 11, 15),  # right biased
    (2, 4, 8, 12),   # compressed-left
  ]

  # Seed with one unit interval to maximize number of KT rounds under CAP.
  T = [(0, 1)]

  # Growth model per KT round: new_size = 4 * size + 4
  def round_next_size(sz):
    return 4 * sz + 4

  def max_rounds_within_cap(initial_size, max_rounds):
    sz = initial_size
    done = 0
    for _ in range(max(0, int(max_rounds))):
      nxt = round_next_size(sz)
      if nxt > CAP:
        break
      sz = nxt
      done += 1
    return done, sz

  depth, _ = max_rounds_within_cap(len(T), rounds)

  def _span_info(current_T):
    lo = min(l for l, r in current_T)
    hi = max(r for l, r in current_T)
    delta = hi - lo
    if delta <= 0:
      delta = 1
    return lo, hi, delta

  def apply_round(current_T, starts, do_interleave=False, round_idx=0):
    lo, hi, delta = _span_info(current_T)

    # Build blocks, optionally reversing source order for parity
    blocks = []
    for idx, s in enumerate(starts):
      base = s * delta - lo
      src = current_T[::-1] if (reverse_block_parity and (idx % 2 == 1)) else current_T
      block = [(l + base, r + base) for (l, r) in src]
      blocks.append(block)

    # Merge blocks: interleave or concatenate
    S = []
    if do_interleave and interleave_blocks:
      maxlen = max((len(b) for b in blocks), default=0)
      order = list(range(len(blocks)))
      # occasionally reverse overall block order to break structure
      if round_idx % 2 == 1:
        order.reverse()
      for i in range(maxlen):
        for idx in order:
          blk = blocks[idx]
          if i < len(blk):
            S.append(blk[i])
    else:
      for blk in blocks:
        S.extend(blk)

    # Classic connectors (KT-style), plus an optional longer-range cross to enhance pressure
    s0, s1, s2, s3 = (starts[0], starts[1], starts[2], starts[3]) if len(starts) >= 4 else (2,6,10,14)
    connectors = [
      ((s0 - 1) * delta, (s1 - 1) * delta),            # left cap
      ((s2 + 2) * delta, (s3 + 2) * delta),            # right cap
      ((s0 + 2) * delta, (s2 - 1) * delta),            # cross 1
      ((s1 + 2) * delta, (s3 - 1) * delta),            # cross 2
    ]
    # Add a slightly longer connector every other round to couple distant blocks
    if round_idx % 2 == 0:
      connectors.append(((s0 + 3) * delta, (s3 + 3) * delta))
    # Ensure valid connectors (r > l)
    for a, b in connectors:
      if b > a:
        S.append((a, b))

    return S

  # Stage 1: build spine rounds
  for ridx in range(depth):
    starts = template_bank[ridx % len(template_bank)] if rotate_starts else template_bank[0]
    do_inter = bool(interleave_blocks and (ridx % 2 == 0))
    T = apply_round(T, starts, do_interleave=do_inter, round_idx=ridx)
    if len(T) >= CAP:
      T = T[:CAP]
      break

  # If already near capacity, return baseline
  if len(T) >= CAP - 16:
    # normalize and convert below
    pass
  else:
    # Stage 2: micro-phase -- up to two small delta-scaled rounds to increase FF usage
    def thin_seed(current_T, target_sz):
      n = len(current_T)
      if n == 0 or target_sz <= 0:
        return []
      step = max(1, n // target_sz)
      U = current_T[::step][:target_sz]
      return U

    def build_micro_delta_round(current_T, budget, iter_id=0, alt=False):
      if not current_T or budget <= 6:
        return []

      glo = min(l for l, r in current_T)
      ghi = max(r for l, r in current_T)
      G = max(1, ghi - glo)

      # choose seed size small but dependent on budget and overall size
      seed_target = max(8, min(64, budget // 10, len(current_T) // 200 if len(current_T) >= 200 else 8))
      U = thin_seed(current_T, seed_target)
      if not U:
        return []

      ulo = min(l for l, r in U)

      # construct four fractional windows; apply small iter-dependent shift to break symmetry
      base_shifts = [0.00, 0.02, -0.02, 0.04] if not alt else [0.01, -0.01, 0.03, -0.03]
      window_fracs = [
        (0.12 + base_shifts[0], 0.22 + base_shifts[0]),
        (0.35 + base_shifts[1], 0.45 + base_shifts[1]),
        (0.58 + base_shifts[2], 0.68 + base_shifts[2]),
        (0.80 + base_shifts[3], 0.90 + base_shifts[3]),
      ]
      # clamp fracs into (0,1)
      window_fracs = [(
        max(0.02, min(0.88, a)),
        max(0.10, min(0.98, b))
      ) for (a, b) in window_fracs]

      blocks = []
      for wi, (fa, fb) in enumerate(window_fracs):
        win_lo = glo + int(round(fa * G))
        base = win_lo - ulo
        block = [(l + base, r + base) for (l, r) in U]
        # alternate internal order to diversify overlaps
        if (wi + iter_id + (1 if alt else 0)) % 2 == 0:
          block = list(reversed(block))
        blocks.append(block)

      # Interleave micro-blocks; alternate ordering for variety
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

      # Add deterministic connectors at reduced scale to couple windows
      # Use half-scale delta to keep these local and safe
      delta2 = max(1, G // 2)
      s0, s1, s2, s3 = (2, 6, 10, 14)
      micro_connectors = [
        ((s0 - 1) * delta2 + glo, (s1 - 1) * delta2 + glo),
        ((s2 + 2) * delta2 + glo, (s3 + 2) * delta2 + glo),
        ((s0 + 2) * delta2 + glo, (s2 - 1) * delta2 + glo),
        ((s1 + 2) * delta2 + glo, (s3 - 1) * delta2 + glo),
      ]
      if alt:
        # add a slightly longer connector in alternate micro-round
        micro_connectors.append(((s0 + 3) * delta2 + glo, (s3 + 3) * delta2 + glo))

      for a, b in micro_connectors:
        if b > a:
          micro.append((a, b))

      # add a few sparse caps to raise FF pressure while keeping cliques bounded
      mid = glo + G // 2
      cap_candidates = [
        (glo + int(0.08 * G), glo + int(0.30 * G)),
        (glo + int(0.60 * G), glo + int(0.92 * G)),
        (mid - max(1, G // 32), mid + max(1, G // 32)),
      ]
      for cap in cap_candidates:
        if cap[1] > cap[0]:
          micro.append(cap)

      # enforce budget
      if len(micro) > budget:
        micro = micro[:budget]
      return micro

    # run up to two micro-steps if space allows
    steps = min(max(0, int(phase2_iters)), 2)
    for iter_id in range(steps):
      room = CAP - len(T)
      if room <= 6:
        break
      micro = build_micro_delta_round(T, room, iter_id=iter_id, alt=(iter_id % 2 == 1))
      if not micro:
        break
      avail = CAP - len(T)
      if len(micro) > avail:
        micro = micro[:avail]
      T.extend(micro)

  # Final normalization and integerization
  if not T:
    return []

  min_l = min(l for l, r in T)
  if min_l < 0:
    T = [(l - min_l, r - min_l) for (l, r) in T]

  intervals = []
  for (l, r) in T:
    li = int(round(l))
    ri = int(round(r))
    if ri <= li:
      ri = li + 1
    intervals.append((li, ri))

  # Trim to CAP just in case
  if len(intervals) > CAP:
    intervals = intervals[:CAP]

  return intervals

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()