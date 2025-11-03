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
    # Round-robin rotate among four strong templates to couple colors across scales
    starts = template_bank[ridx % 4]
    # Parity-based interleaving (on even rounds) to enhance FF mixing without large omega increase
    do_interleave = (ridx % 2 == 0)
    T = apply_round(T, starts, do_interleave=do_interleave)

    # On the final spine round, add a single long-range cross connector (cross4) to couple colors across distant spans
    if ridx == depth - 1 and len(T) <= CAP - 1:
      # compute span based on the just-built T (safe, deterministic)
      lo_final = min(l for l, r in T)
      hi_final = max(r for l, r in T)
      delta_final = max(1, hi_final - lo_final)
      s0, s3 = starts[0], starts[3]
      T.append(((s0 + 4) * delta_final, (s3 + 4) * delta_final))

  # If we already reached or are near the cap, skip phase 2 safely.
  if len(T) >= CAP - 16:
    return T

  # Inject a small set of long-range tail caps (capacity-safe) to increase FF pressure
  lo_span = min(l for l, r in T)
  hi_span = max(r for l, r in T)
  span = max(1, hi_span - lo_span)
  def cap_at(a_frac, b_frac):
    L = lo_span + max(1, int(round(a_frac * span)))
    R = lo_span + max(1, int(round(b_frac * span)))
    if R <= L:
      R = L + 1
    return (L, R)
  tail_caps = [cap_at(0.08, 0.60), cap_at(0.25, 0.75), cap_at(0.75, 0.92)]
  room_caps = CAP - len(T)
  if room_caps > 0:
    for c in tail_caps:
      if len(T) >= CAP:
        break
      T.append(c)

  # Stage 2: delta2-driven micro round using a thin seed of T to raise FF without
  # pushing omega too high. Deterministic, capacity-bounded, and single-pass by default.
  def build_micro_delta_round(current_T, budget):
    if not current_T or budget <= 8:
      return []

    # Global span
    glo = min(l for l, r in current_T)
    ghi = max(r for l, r in current_T)
    G = max(1, ghi - glo)

    # Sample a thin seed from current_T, evenly spaced to diversify structure.
    # Keep it small and bounded by budget.
    max_seed_by_budget = max(8, min(32, (budget - 4) // 4))
    max_seed_by_len = max(8, len(current_T) // 300)
    seed_sz = min(max_seed_by_budget, max_seed_by_len, 32)
    seed_sz = max(8, seed_sz)
    stride = max(1, len(current_T) // seed_sz)
    U = [current_T[i] for i in range(0, len(current_T), stride)][:seed_sz]
    if not U:
      return []

    ulo = min(l for l, r in U)

    # Four windows across the global span (fractions chosen to stay well inside [glo, ghi])
    # This emulates a small KT-like four-block round at a reduced scale.
    window_fracs = [(0.12, 0.22), (0.35, 0.45), (0.58, 0.68), (0.80, 0.90)]

    # Build translated micro-blocks aligned to these windows
    blocks = []
    for (fa, fb) in window_fracs:
      win_lo = glo + int(round(fa * G))
      win_hi = glo + int(round(fb * G))
      base = win_lo - ulo
      block = [(l + base, r + base) for (l, r) in U]
      blocks.append(block)

    # Interleave micro-blocks to maximize FF mixing
    micro = []
    maxlen = max(len(b) for b in blocks)
    for i in range(maxlen):
      for blk in blocks:
        if i < len(blk):
          micro.append(blk[i])

    # Deterministic connectors across windows (fractional-span analog of KT caps)
    connectors = [
      (glo + int(round(0.08 * G)), glo + int(round(0.30 * G))),  # left cap
      (glo + int(round(0.60 * G)), glo + int(round(0.92 * G))),  # right cap
      (glo + int(round(0.26 * G)), glo + int(round(0.56 * G))),  # cross 1
      (glo + int(round(0.44 * G)), glo + int(round(0.78 * G))),  # cross 2
    ]
    micro.extend(connectors)

    # Trim to available budget
    if len(micro) > budget:
      micro = micro[:budget]
    return micro

  # Respect requested phase2_iters but guard to at most two micro-rounds by default.
  steps = min(max(0, int(phase2_iters)), 2)
  for _ in range(steps):
    room = CAP - len(T)
    if room <= 8:
      break
    micro = build_micro_delta_round(T, room)
    if not micro:
      break
    # Capacity guard and append
    avail = CAP - len(T)
    if len(micro) > avail:
      micro = micro[:avail]
    T.extend(micro)

  # Optional second micro-phase with an alternate window family to diversify coupling
  def build_micro_delta_round_alt(current_T, budget):
    if not current_T or budget <= 8:
      return []
    glo = min(l for l, r in current_T)
    ghi = max(r for l, r in current_T)
    G = max(1, ghi - glo)

    # Thin seed: slightly smaller to control omega at alt scale
    seed_sz = max(8, min(28, len(current_T) // 320))
    stride = max(1, len(current_T) // seed_sz)
    U = [current_T[i] for i in range(0, len(current_T), stride)][:seed_sz]
    if not U:
      return []
    ulo = min(l for l, r in U)

    # Alternate windows (distinct placement)
    window_fracs = [(0.05, 0.15), (0.28, 0.38), (0.60, 0.70), (0.82, 0.92)]
    blocks = []
    for (fa, fb) in window_fracs:
      win_lo = glo + int(round(fa * G))
      base = win_lo - ulo
      blk = [(l + base, r + base) for (l, r) in U]
      blocks.append(blk)

    micro = []
    maxlen = max(len(b) for b in blocks)
    order = [0, 2, 1, 3]  # fixed interleave order to break symmetry
    for i in range(maxlen):
      for idx in order:
        b = blocks[idx]
        if i < len(b):
          micro.append(b[i])

    # Fractional connectors tuned for alt windows
    connectors = [
      (glo + int(round(0.10 * G)), glo + int(round(0.28 * G))),
      (glo + int(round(0.52 * G)), glo + int(round(0.88 * G))),
      (glo + int(round(0.22 * G)), glo + int(round(0.62 * G))),
      (glo + int(round(0.40 * G)), glo + int(round(0.74 * G))),
    ]
    micro.extend(connectors)

    if len(micro) > budget:
      micro = micro[:budget]
    return micro

  room2 = CAP - len(T)
  if room2 > 8:
    microB = build_micro_delta_round_alt(T, room2)
    if microB:
      avail2 = CAP - len(T)
      if len(microB) > avail2:
        microB = microB[:avail2]
      T.extend(microB)

  return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()