# EVOLVE-BLOCK-START

def construct_intervals(rounds=6,
                        rotate_starts=True,
                        reverse_block_parity=True,
                        interleave_blocks=True,
                        phase2_iters=2,
                        seed_count=1):
  """
  Hex-spine with deterministic rotation, enforced interleaving, dual CAP-aware micro phases,
  and bounded long-range connectors. Returns intervals (l, r) in FirstFit presentation order.
  """

  # Hard capacity guard to keep total intervals < 10000
  CAP = 9800

  # Six-template start bank (hex spine)
  # The two extra templates beyond the classic four diversify color coupling across rounds.
  template_bank6 = [
    (2, 6, 10, 14),  # T1: classic KT
    (1, 5, 9, 13),   # T2: left-shifted
    (3, 7, 11, 15),  # T3: right-shifted
    (4, 8, 12, 16),  # T4: stretched-right
    (2, 4, 8, 12),   # T5: compressed-left pair
    (5, 9, 13, 17),  # T6: stretched-far right
  ]

  # Deterministic seed for template rotation (no randomness; reproducible)
  BASE_SEED = 131071

  # Seed initialization: one unit interval; multi-seed allowed but constrained
  if seed_count <= 1:
    T = [(0, 1)]
  else:
    seeds = min(4, max(1, int(seed_count)))
    step = 3
    T = [(i * step, i * step + 1) for i in range(seeds)]

  # Helpers
  def _span_delta(current_T):
    lo = min(l for l, r in current_T)
    hi = max(r for l, r in current_T)
    delta = hi - lo
    if delta <= 0:
      delta = 1
    return lo, hi, int(delta)

  # Classic KT connectors (four caps/crosses), scale-preserving
  def _append_connectors(S, starts, delta):
    s0, s1, s2, s3 = starts
    S.append(((s0 - 1) * delta, (s1 - 1) * delta))  # left cap
    S.append(((s2 + 2) * delta, (s3 + 2) * delta))  # right cap
    S.append(((s0 + 2) * delta, (s2 - 1) * delta))  # cross 1
    S.append(((s1 + 2) * delta, (s3 - 1) * delta))  # cross 2)

  # Build a spine round: always interleave blocks (mandated), optional reverse of block order
  def _apply_spine_round(current_T, starts, ridx, enforce_interleave=True, allow_reverse=True):
    lo, hi, delta = _span_delta(current_T)

    # Build four translated copies
    blocks = []
    for s in starts:
      base = s * delta - lo
      blk = [(l + base, r + base) for (l, r) in current_T]
      blocks.append(blk)

    # Assemble S with strong interleaving to maximize color mixing
    S = []
    order = list(range(4))
    if allow_reverse and (ridx % 2 == 1):
      order.reverse()
    if enforce_interleave:
      maxlen = max(len(b) for b in blocks)
      for i in range(maxlen):
        for idx in order:
          blk = blocks[idx]
          if i < len(blk):
            S.append(blk[i])
    else:
      for idx in order:
        S.extend(blocks[idx])

    # Classic connectors
    _append_connectors(S, starts, delta)
    return S

  # Tiny per-round density booster: add a handful of short pins near edges, CAP-guarded
  def _edge_pins(current_T, max_add=6):
    if max_add <= 0 or not current_T:
      return []
    lo, hi, G = _span_delta(current_T)
    if G <= 4:
      return []
    # Place pins in low-density edge windows to avoid inflating omega
    # Fractions chosen away from core spine windows
    fracs = [0.03, 0.07, 0.93, 0.97]
    eps = max(1, G // 512)
    pins = []
    for f in fracs:
      a = lo + int(round(f * G))
      pins.append((a, a + eps))
    return pins[:max_add]

  # Predictive spine growth: sz -> 4*sz + 4 per round
  def _max_rounds_within_cap(initial_size, max_rounds):
    sz = initial_size
    done = 0
    for k in range(max(0, int(max_rounds))):
      nxt = 4 * sz + 4  # KT recurrence
      if nxt > CAP:
        break
      sz = nxt
      done += 1
    return done, sz

  # Stage 1: Hex spine with deterministic rotation and mandated interleaving
  depth, _ = _max_rounds_within_cap(len(T), rounds)
  for ridx in range(depth):
    # Deterministic template selection (six-template rotation)
    if rotate_starts:
      h = (BASE_SEED * (ridx + 1) + 977 * ridx) % 6
      starts = template_bank6[h]
    else:
      starts = template_bank6[0]
    # Enforce interleaving irrespective of the incoming flag to maximize pressure
    T = _apply_spine_round(T, starts, ridx,
                           enforce_interleave=True,
                           allow_reverse=True if reverse_block_parity else False)
    # Lightweight edge pins to slightly densify without raising omega (strictly bounded)
    if len(T) < CAP - 16:
      pins = _edge_pins(T, max_add=4 if ridx % 2 == 0 else 2)
      room = CAP - len(T)
      if room > 0 and pins:
        T.extend(pins[:room])

    if len(T) >= CAP:
      return T[:CAP]

  # Early return if nearly full
  if len(T) >= CAP - 16:
    return T[:CAP]

  # Stage 1.5: Long-range deterministic connectors (guarded)
  def _global_connectors(current_T, limit=6):
    lo, hi, G = _span_delta(current_T)
    con_fracs = [(0.05, 0.37), (0.18, 0.62), (0.26, 0.78), (0.44, 0.86), (0.60, 0.92), (0.08, 0.30)]
    out = []
    for (fa, fb) in con_fracs[:limit]:
      a = lo + int(round(fa * G))
      b = lo + int(round(fb * G))
      if b > a:
        out.append((a, b))
    return out

  # Reserve a fraction of remaining CAP for connectors and micro phases
  remain = CAP - len(T)
  if remain > 0:
    # Connector budget clamped small to avoid omega inflation
    conn_budget = min(10, max(4, remain // 1000))
    conns = _global_connectors(T, limit=conn_budget)
    add = conns[:remain]
    T.extend(add)
  if len(T) >= CAP - 8:
    return T[:CAP]

  # Stage 2: Dual micro phases with explicit CAP budgeting
  # Allocate a fraction of the remaining capacity to micro passes
  remain = CAP - len(T)
  if remain <= 8:
    return T[:CAP]

  micro_total = remain
  # Reserve budgets: A gets 55%, B gets 45% (rounded)
  budget_A = int(0.55 * micro_total)
  budget_B = micro_total - budget_A

  # Thin seed
  def _thin_seed(current_T, max_seed):
    n = len(current_T)
    if n == 0 or max_seed <= 0:
      return []
    step = max(1, n // max_seed)
    return [current_T[i] for i in range(0, n, step)][:max_seed]

  def _build_micro(current_T, budget, window_fracs, connectors=True, reverse=False):
    if not current_T or budget <= 0:
      return []
    glo = min(l for l, r in current_T)
    ghi = max(r for l, r in current_T)
    G = max(1, ghi - glo)

    # Seed size: small to protect omega, yet provide mixing
    seed_sz = max(12, min(40, len(current_T) // 250))
    U = _thin_seed(current_T, seed_sz)
    if not U:
      return []
    ulo = min(l for l, r in U)

    # Build translated micro-blocks aligned to windows
    blocks = []
    for (fa, fb) in window_fracs:
      win_lo = glo + int(round(fa * G))
      base = win_lo - ulo
      blk = [(l + base, r + base) for (l, r) in U]
      blocks.append(blk)
    if reverse:
      blocks = list(reversed(blocks))

    # Interleave micro blocks fully
    micro = []
    maxlen = max(len(b) for b in blocks)
    for i in range(maxlen):
      for blk in blocks:
        if i < len(blk):
          micro.append(blk[i])

    # Add fractional connectors across the micro windows (short)
    if connectors:
      micro_connectors = [
        (glo + int(round(0.08 * G)), glo + int(round(0.30 * G))),
        (glo + int(round(0.26 * G)), glo + int(round(0.56 * G))),
        (glo + int(round(0.44 * G)), glo + int(round(0.78 * G))),
        (glo + int(round(0.60 * G)), glo + int(round(0.92 * G))),
      ]
      for a, b in micro_connectors:
        if b > a:
          micro.append((a, b))

    # Trim to budget
    if len(micro) > budget:
      micro = micro[:budget]
    return micro

  # Micro-pass A: left/inner windows (denser toward the core)
  windows_A = [(0.12, 0.22), (0.36, 0.46), (0.58, 0.68), (0.80, 0.90)]
  microA = _build_micro(T, budget_A, windows_A, connectors=True, reverse=False)
  if microA:
    room = CAP - len(T)
    if len(microA) > room:
      microA = microA[:room]
    T.extend(microA)

  # Micro-pass B: right/outer windows (shifted to avoid cliques and overlap cores)
  remain = CAP - len(T)
  if remain > 0:
    budget_B = min(budget_B, remain)
  windows_B = [(0.10, 0.20), (0.30, 0.40), (0.62, 0.72), (0.84, 0.94)]
  microB = _build_micro(T, budget_B, windows_B, connectors=True, reverse=True)
  if microB:
    room = CAP - len(T)
    if len(microB) > room:
      microB = microB[:room]
    T.extend(microB)

  # Final tiny pin fill if space remains (keep very small to avoid omega growth)
  remain = CAP - len(T)
  if remain > 0:
    lo, hi, G = _span_delta(T)
    eps = max(1, G // 1024)
    pins = []
    # Deterministic positions away from core windows
    for f in (0.06, 0.24, 0.76, 0.94):
      a = lo + int(round(f * G))
      pins.append((a, a + eps))
    if pins:
      T.extend(pins[:remain])

  # Final capacity trim
  if len(T) > CAP:
    T = T[:CAP]
  return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()