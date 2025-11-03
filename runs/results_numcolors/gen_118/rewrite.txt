# EVOLVE-BLOCK-START

def construct_intervals(rounds=6,
                        rotate_starts=True,
                        reverse_block_parity=True,
                        interleave_blocks=True,
                        phase2_iters=0):
  """
  Enhanced deterministic KT-style spine and capped micro-extensions.

  Signature preserved for compatibility with previous evaluator.
  Returns a list of open intervals (l, r) in the order presented to FirstFit.
  """

  # Hard capacity guard to keep the total count < 10000
  CAP = 9800

  # Expanded template bank (8 templates) of 4-start patterns.
  template_bank = [
    (2, 6, 10, 14),  # classic
    (1, 5, 9, 13),   # left-shift
    (3, 7, 11, 15),  # right-shift
    (2, 4, 8, 12),   # compressed left pair
    (3, 5, 9, 13),   # gentle left pack
    (1, 7, 11, 15),  # wide skew
    (2, 8, 10, 12),  # inner symmetric
    (4, 6, 8, 10),   # tight middle
  ]

  # Initial seed (single unit interval as before to permit depth)
  T = [(0.0, 1.0)]

  # Growth rule: each KT-style round roughly maps size -> 4*size + 4.
  def next_round_size(sz):
    return 4 * sz + 4

  # Compute maximal feasible rounds under CAP
  def max_rounds_allowable(initial_size, max_rounds):
    sz = initial_size
    used = 0
    for _ in range(int(max_rounds)):
      nxt = next_round_size(sz)
      if nxt > CAP:
        break
      sz = nxt
      used += 1
    return used, sz

  depth, _ = max_rounds_allowable(len(T), rounds)

  # Thin deterministic sampler: uniformly subsample at most max_seed items
  def thin_seed(current_T, max_seed):
    n = len(current_T)
    if n == 0 or max_seed <= 0:
      return []
    step = max(1, n // max_seed)
    return current_T[::step][:max_seed]

  # Apply one KT-style spine round with given start pattern and interleaving mode.
  # We also insert a small number of "spine blockers" per block to occupy early FF colors.
  def apply_spine_round(current_T, starts, round_idx, do_interleave):
    lo = min(l for l, r in current_T)
    hi = max(r for l, r in current_T)
    delta = hi - lo
    if delta <= 0:
      delta = 1.0

    # Build translated blocks
    blocks = []
    for s in starts:
      base = s * delta - lo
      block = [(l + base, r + base) for (l, r) in current_T]
      blocks.append(block)

    # Insert tiny in-block spine blockers (longer intervals placed first in each block)
    # They are local to each translated block (slightly extend the block span).
    # Keep count small to avoid raising the global clique number too much.
    spine_blockers = []
    blocker_len = max(1.0, delta * 0.6)  # somewhat long compared to block internal width
    for idx, s in enumerate(starts):
      # place blocker near the block's left edge (deterministic offset)
      base = s * delta - lo
      b_lo = base - 0.05 * delta + 0.02 * idx
      b_hi = b_lo + blocker_len
      spine_blockers.append((b_lo, b_hi))

    # Produce final sequence S for this round with the chosen interleaving policy.
    S = []
    if do_interleave:
      # two-round adaptive parity: rounds in {0,1} forward interleave, {2,3} reverse, cycle
      cadence = round_idx % 4
      maxlen = max(len(b) for b in blocks)
      if cadence in (0, 1):
        # forward interleave
        for i in range(maxlen):
          for blk in blocks:
            if i < len(blk):
              S.append(blk[i])
      else:
        # reverse interleave
        for i in range(maxlen):
          for blk in reversed(blocks):
            if i < len(blk):
              S.append(blk[i])
    else:
      # sequential block order
      for blk in blocks:
        S.extend(blk)

    # Prepend spine blockers so they appear before the block intervals in each round,
    # pressuring FirstFit to use small colors early.
    S = spine_blockers + S

    # Classic KT-style connectors adapted to the current delta and starts.
    s0, s1, s2, s3 = starts
    connectors = [
      ((s0 - 1) * delta, (s1 - 1) * delta),  # left cap
      ((s2 + 2) * delta, (s3 + 2) * delta),  # right cap
      ((s0 + 2) * delta, (s2 - 1) * delta),  # cross 1
      ((s1 + 2) * delta, (s3 - 1) * delta),  # cross 2
    ]
    S.extend(connectors)

    return S

  # Run spine rounds cycling through template bank and interleaving cadence.
  for ridx in range(depth):
    starts = template_bank[ridx % len(template_bank)]
    # enable interleaving for most rounds to increase mixing; still deterministic
    do_interleave = bool(interleave_blocks)
    T = apply_spine_round(T, starts, ridx, do_interleave=do_interleave)
    # Safety cap: if we near capacity, stop expanding spine
    if len(T) > CAP - 64:
      break

  # If we already reached or are very near the cap, return T
  if len(T) >= CAP - 16:
    return T

  # Stage 2: improved micro-extension phase.
  # We'll attempt up to two delta2 rounds plus one tighter delta2' round,
  # each strictly bounded by remaining capacity.

  def micro_round(current_T, round_id, budget, scale_factor=0.5, extra_long_caps=3):
    """
    scale_factor: fraction of global span used as local delta2 scale.
    extra_long_caps: number of sparse long caps added to overlap many active colors.
    """
    if budget <= 0 or not current_T:
      return []

    glo = min(l for l, r in current_T)
    ghi = max(r for l, r in current_T)
    G = max(1.0, ghi - glo)

    # delta2 at controlled scale (integer-ish where possible)
    delta2 = max(1.0, max(1.0, scale_factor * G))

    # Thin deterministic seed; target a modest seed size to avoid large omega
    per_block_target = max(12, min(96, budget // 10))
    U = thin_seed(current_T, per_block_target)
    if not U:
      return []

    # Build four micro-blocks using the classic 4-start pattern from the template bank
    starts = template_bank[round_id % len(template_bank)]
    ulo = min(l for l, r in U)
    blocks = []
    for s in starts:
      base = s * delta2 + glo - ulo
      block = [(l + base, r + base) for (l, r) in U]
      # deterministic micro reversal to diversify interaction pattern
      if (s % 2) == (round_id % 2):
        block = list(reversed(block))
      blocks.append(block)

    # Interleave blocks using the same two-round cadence (round_id influences order)
    micro = []
    cadence = round_id % 4
    maxlen = max(len(b) for b in blocks)
    if cadence in (0, 1):
      for i in range(maxlen):
        for blk in blocks:
          if i < len(blk):
            micro.append(blk[i])
    else:
      for i in range(maxlen):
        for blk in reversed(blocks):
          if i < len(blk):
            micro.append(blk[i])

    # Deterministic connectors at delta2 scale (localized to glo..ghi region)
    s0, s1, s2, s3 = starts
    connectors = [
      (glo + (s0 - 1) * delta2, glo + (s1 - 1) * delta2),
      (glo + (s2 + 2) * delta2, glo + (s3 + 2) * delta2),
      (glo + (s0 + 2) * delta2, glo + (s2 - 1) * delta2),
      (glo + (s1 + 2) * delta2, glo + (s3 - 1) * delta2),
    ]
    micro.extend(connectors)

    # Add a few sparse long caps that overlap many blocks — tuned to be 'safe' yet pressuring FF.
    caps = []
    for k in range(extra_long_caps):
      frac_lo = 0.05 + 0.2 * k
      frac_hi = min(0.95, frac_lo + 0.5)
      cap_lo = glo + frac_lo * G
      cap_hi = glo + frac_hi * G
      if cap_hi > cap_lo:
        caps.append((cap_lo, cap_hi))
    micro.extend(caps)

    # Trim micro to budget deterministically
    if len(micro) > budget:
      micro = micro[:budget]
    return micro

  # Two delta2 micro rounds
  remaining = CAP - len(T)
  # allocate roughly half of remaining to first micro, rest to second
  if remaining > 16:
    b1 = max(8, remaining // 2)
    mr1 = micro_round(T, round_id=0, budget=b1, scale_factor=0.45, extra_long_caps=3)
    if mr1:
      room = CAP - len(T)
      if len(mr1) > room:
        mr1 = mr1[:room]
      T.extend(mr1)

  remaining = CAP - len(T)
  if remaining > 16:
    b2 = max(8, remaining)
    mr2 = micro_round(T, round_id=1, budget=b2, scale_factor=0.35, extra_long_caps=3)
    if mr2:
      room = CAP - len(T)
      if len(mr2) > room:
        mr2 = mr2[:room]
      T.extend(mr2)

  # Optional tight micro (delta2') to push FF late: use a smaller scale and tiny seed.
  remaining = CAP - len(T)
  if remaining > 32:
    b3 = min(remaining, 256)
    # very fine-scale micro extension: smaller seed, tiny scale to avoid big span growth
    mr3 = micro_round(T, round_id=2, budget=b3, scale_factor=0.12, extra_long_caps=2)
    if mr3:
      room = CAP - len(T)
      if len(mr3) > room:
        mr3 = mr3[:room]
      T.extend(mr3)

  # Final safety trim (if anything overshot due to arithmetic)
  if len(T) > CAP:
    T = T[:CAP]

  return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()