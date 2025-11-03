# EVOLVE-BLOCK-START

def construct_intervals(rounds=6,
                        rotate_starts=True,
                        reverse_block_parity=True,
                        interleave_blocks=True,
                        phase2_iters=0):
  """
  Deterministic KT-style spine with a safeguarded two-phase scaffold.

  Parameters:
    rounds (int): main expansion depth; near 6 yields ~9556 intervals for a single-seed KT spine.
    rotate_starts (bool): available but disabled in the safeguarded spine mode.
    reverse_block_parity (bool): available but disabled in the safeguarded spine mode.
    interleave_blocks (bool): available but disabled in the safeguarded spine mode.
    phase2_iters (int): requested micro iterations; actual application is strictly capacity- and safety-gated.

  Returns:
    intervals: list of (l, r) integer tuples, open intervals, in FF presentation order.
  """

  # Capacity guard relaxed to permit a thin micro-phase while staying < 10000
  CAP = 9800

  # Stable KT-spine start pattern (proven strong in prior runs)
  spine_starts = [2, 6, 10, 14]

  # Optional template bank for exploratory use in micro/secondary phase
  # (kept dormant by default to preserve the strong baseline).
  template_bank = [
    [2, 6, 10, 14],  # A: classic KT
    [1, 5, 9, 13],   # B: left-shifted
    [3, 7, 11, 15],  # C: right-shifted
    [2, 4, 8, 12],   # D: compressed left pair
    [3, 5, 9, 13],   # E: gentle left pack
    [1, 7, 11, 15],  # F: wide skew
    [2, 8, 10, 12],  # G: inner symmetric
    [4, 6, 8, 10],   # H: tight middle
  ]

  # Seed with one unit interval to allow six KT rounds within CAP.
  T = [(0, 1)]

  # Predictive size accounting to cap the number of full rounds.
  # KT growth per round: size -> 4*size + 4
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

  # Stage 1: KT spine, un-interleaved, classic connectors, capacity-safe depth.
  depth, size_est = max_rounds_within_cap(len(T), rounds)

  def apply_round(current_T, starts, do_interleave=False):
    # Compute span and delta
    lo = min(l for l, r in current_T)
    hi = max(r for l, r in current_T)
    delta = hi - lo
    if delta <= 0:
      delta = 1

    # Build four translated blocks; optionally interleave to enhance color mixing
    blocks = []
    for s in starts:
      base = s * delta - lo
      block = [(l + base, r + base) for (l, r) in current_T]
      blocks.append(block)

    # Interleave on demand (even rounds), otherwise keep sequential (odd rounds)
    S = []
    if do_interleave:
      maxlen = max(len(b) for b in blocks)
      for i in range(maxlen):
        for blk in blocks:
          if i < len(blk):
            S.append(blk[i])
    else:
      for blk in blocks:
        S.extend(blk)

    # Classic connectors that preserve the desired FF pressure without
    # blowing up omega when used with spine_starts:
    s0, s1, s2, s3 = starts
    connectors = [
      ((s0 - 1) * delta, (s1 - 1) * delta),
      ((s2 + 2) * delta, (s3 + 2) * delta),
      ((s0 + 2) * delta, (s2 - 1) * delta),
      ((s1 + 2) * delta, (s3 - 1) * delta),
    ]
    S.extend(connectors)
    return S

  for ridx in range(depth):
    # Alternate interleaving on even rounds to increase FF pressure without large omega spikes
    T = apply_round(T, spine_starts, do_interleave=(ridx % 2 == 0))

  # If we already reached or are near the cap, skip phase 2 safely.
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

    # Use a round-dependent reduced scale to avoid expanding span while adding coupling.
    # round_id: 0 -> G//2, 1 -> G//3, >=2 -> G//4
    div = 2 if round_id == 0 else (3 if round_id == 1 else 4)
    delta2 = max(1, G // div)

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
    starts = spine_starts
    blocks = []
    ulo = min(l for l, r in U)
    for s in starts:
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
    s0, s1, s2, s3 = starts
    connectors = [
      ((s0 - 1) * delta2 + glo, (s1 - 1) * delta2 + glo),  # left cap (localized)
      ((s2 + 2) * delta2 + glo, (s3 + 2) * delta2 + glo),  # right cap
      ((s0 + 2) * delta2 + glo, (s2 - 1) * delta2 + glo),  # cross 1
      ((s1 + 2) * delta2 + glo, (s3 - 1) * delta2 + glo),  # cross 2
      ((s0 + 3) * delta2 + glo, (s3 + 3) * delta2 + glo),  # cross3 (longer range)
      ((s0 + 4) * delta2 + glo, (s3 + 4) * delta2 + glo),  # cross4 (additional long-range coupler)
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

  # Optional third micro round at an even smaller scale to add gentle late pressure
  remaining = CAP - len(T)
  mr3 = micro_round(T, round_id=2, budget=max(0, remaining))
  if mr3:
    room = CAP - len(T)
    if len(mr3) > room:
      mr3 = mr3[:room]
    T.extend(mr3)

  return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()