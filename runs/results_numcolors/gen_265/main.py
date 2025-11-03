# EVOLVE-BLOCK-START

def construct_intervals(rounds=6,
                        rotate_starts=True,
                        reverse_block_parity=True,
                        interleave_blocks=True,
                        phase2_iters=2):
  """
  Rotating four-template KT spine with an enabled micro-phase to boost FirstFit.

  Parameters:
    rounds (int): main expansion depth; near 6 yields ~9556 intervals for a single-seed KT spine.
    rotate_starts (bool): available but disabled in the safeguarded spine mode.
    reverse_block_parity (bool): available but disabled in the safeguarded spine mode.
    interleave_blocks (bool): available but disabled in the safeguarded spine mode.
    phase2_iters (int): requested micro iterations; actual application is strictly capacity- and safety-gated.

  Returns:
    intervals: list of (l, r) integer tuples, open intervals, in FF presentation order.
  """

  # Hard capacity guard raised to allow micro-phase after six KT rounds
  CAP = 9800

  # Stable KT-spine start pattern (proven strong in prior runs)
  spine_starts = [2, 6, 10, 14]

  # Rotating template bank for the main spine (round-robin over the first four)
  template_bank = [
    [2, 6, 10, 14],  # T1: classic KT
    [1, 5, 9, 13],   # T2: left-shifted
    [3, 7, 11, 15],  # T3: right-shifted
    [4, 8, 12, 16],  # T4: stretched right (new)
    [2, 4, 8, 12],   # aux: compressed left pair
    [3, 5, 9, 13],   # aux: gentle left pack
    [1, 7, 11, 15],  # aux: wide skew
    [2, 8, 10, 12],  # aux: inner symmetric
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
    # Round-robin rotate among four strong templates to couple colors across scales.
    # Enable parity-based interleaving: even rounds interleave, odd rounds sequential.
    starts = template_bank[ridx % 4]
    do_interleave = (ridx % 2 == 0)
    T = apply_round(T, starts, do_interleave=do_interleave)

  # If we already reached or are near the cap, skip phase 2 safely.
  if len(T) >= CAP - 16:
    return T

  # Cross-scale connector pack: a small set of long but sparse connectors placed
  # at fractional positions of the current span. This couples many active colors
  # while keeping omega modest due to sparse use.
  glo = min(l for l, r in T)
  ghi = max(r for l, r in T)
  G = max(1, ghi - glo)
  connector_pack = [
    (glo + int(0.05 * G), glo + int(0.35 * G)),
    (glo + int(0.30 * G), glo + int(0.70 * G)),
    (glo + int(0.65 * G), glo + int(0.95 * G)),
    (glo + int(0.18 * G), glo + int(0.84 * G)),
  ]
  # Capacity-guarded append
  for a, b in connector_pack:
    if len(T) >= CAP - 8:
      break
    if b > a:
      T.append((a, b))

  # Stage 2: delta2-driven micro round using a thin seed of T to raise FF without
  # pushing omega too high. Deterministic, capacity-bounded, and single-pass by default.
  def build_micro_delta_round(current_T, budget, iter_id=0):
    if not current_T or budget <= 8:
      return []

    # Global span
    glo = min(l for l, r in current_T)
    ghi = max(r for l, r in current_T)
    G = max(1, ghi - glo)

    # Sample a thin seed from current_T, evenly spaced to diversify structure.
    # Slightly higher ceiling than before to enhance coupling.
    max_seed_by_budget = max(8, min(40, (budget - 4) // 4))
    max_seed_by_len = max(8, len(current_T) // 280)
    seed_sz = min(max_seed_by_budget, max_seed_by_len, 40)
    seed_sz = max(8, seed_sz)
    stride = max(1, len(current_T) // seed_sz)
    U = [current_T[i] for i in range(0, len(current_T), stride)][:seed_sz]
    if not U:
      return []

    ulo = min(l for l, r in U)

    # Two alternating window families: A (inner windows) and B (shifted, slightly wider).
    if iter_id % 2 == 0:
      window_fracs = [(0.12, 0.22), (0.35, 0.45), (0.58, 0.68), (0.80, 0.90)]
    else:
      window_fracs = [(0.08, 0.18), (0.30, 0.42), (0.55, 0.69), (0.82, 0.94)]

    # Build translated micro-blocks aligned to these windows
    blocks = []
    for (fa, fb) in window_fracs:
      win_lo = glo + int(round(fa * G))
      base = win_lo - ulo
      block = [(l + base, r + base) for (l, r) in U]
      # Break symmetry: reverse every other block on odd iterations
      if iter_id % 2 == 1 and int(round(fa * 100)) // 5 % 2 == 1:
        block = list(reversed(block))
      blocks.append(block)

    # Interleave micro-blocks: forward on even iterations, reverse block order on odd.
    micro = []
    maxlen = max(len(b) for b in blocks)
    block_order = list(range(len(blocks)))
    if iter_id % 2 == 1:
      block_order.reverse()
    for i in range(maxlen):
      for idx in block_order:
        blk = blocks[idx]
        if i < len(blk):
          micro.append(blk[i])

    # Deterministic connectors across windows (fractional-span analog of KT caps)
    connectors = [
      (glo + int(round(0.08 * G)), glo + int(round(0.30 * G))),  # left cap
      (glo + int(round(0.60 * G)), glo + int(round(0.92 * G))),  # right cap
      (glo + int(round(0.26 * G)), glo + int(round(0.56 * G))),  # cross 1
      (glo + int(round(0.44 * G)), glo + int(round(0.78 * G))),  # cross 2
    ]
    # On odd iterations, add a longer-range cross to strengthen cross-scale coupling.
    if iter_id % 2 == 1:
      connectors.append((glo + int(round(0.18 * G)), glo + int(round(0.84 * G))))
    micro.extend(connectors)

    # Trim to available budget
    if len(micro) > budget:
      micro = micro[:budget]
    return micro

  # Respect requested phase2_iters but guard to at most two micro-rounds by default.
  steps = min(max(0, int(phase2_iters)), 2)
  for iter_id in range(steps):
    room = CAP - len(T)
    if room <= 8:
      break
    micro = build_micro_delta_round(T, room, iter_id=iter_id)
    if not micro:
      break
    # Capacity guard and append
    avail = CAP - len(T)
    if len(micro) > avail:
      micro = micro[:avail]
    T.extend(micro)

  return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()