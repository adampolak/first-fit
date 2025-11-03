# EVOLVE-BLOCK-START

def construct_intervals(seed_count=1):
  """
  Construct a sequence of open intervals (l, r) to maximize FirstFit colors
  under an omega (clique number) cap near 10. Deterministic, CAP-guarded.

  Returns:
    intervals: list of (l, r) integer tuples in FirstFit presentation order.
  """

  # -------------------- Global parameters (CAP and budgets) --------------------
  CAP = 9800                  # hard capacity guard
  BASE_SEED = 0x9E3779B97F4A7C15  # deterministic seed for template/parity selection
  SPINE_ROUNDS_MAX = 6        # target spine depth
  SPINE_RESERVE = 240         # reserve for micro phases and connectors (enables 6 rounds)
  SPINE_ALWAYS_INTERLEAVE = True
  LR_CONNECTORS_ENABLED = True
  LR_CONNECTOR_BUDGET = 6     # max long-range connectors total (CAP-guarded)
  SPINE_PIN_BUDGET = 18       # total short pins injected across spine rounds (CAP-guarded)

  # Micro-phase budgets (fractions of remaining CAP after the spine)
  PHASE2_BUDGET_FRACTION_A = 0.60
  PHASE2_BUDGET_FRACTION_B = 0.40

  # -------------------- Template bank (six templates) -------------------------
  # Expanded template bank improves cross-round mixing without raising omega.
  template_bank = [
    (2, 6, 10, 14),  # T1: classic KT
    (1, 5, 9, 13),   # T2: left-shifted
    (3, 7, 11, 15),  # T3: right-shifted
    (4, 8, 12, 16),  # T4: stretched-right
    (2, 4, 8, 12),   # T5: compressed left pair
    (3, 5, 9, 13),   # T6: gentle left pack
  ]

  # -------------------- Helpers --------------------
  def _span_delta(T):
    lo = min(l for l, _ in T)
    hi = max(r for _, r in T)
    delta = hi - lo
    if delta <= 0:
      delta = 1
    return lo, hi, delta

  def _predict_next_size(sz):
    # KT growth upper bound per round: sz -> 4*sz + 4 (four translated blocks + 4 connectors)
    return 4 * sz + 4

  def _mix64(x):
    # Deterministic 64-bit mixer
    x &= (1 << 64) - 1
    x ^= (x >> 30)
    x = (x * 0xBF58476D1CE4E5B9) & ((1 << 64) - 1)
    x ^= (x >> 27)
    x = (x * 0x94D049BB133111EB) & ((1 << 64) - 1)
    x ^= (x >> 31)
    return x

  def _det_choice(seed, ridx, mod):
    return _mix64(seed ^ (ridx * 0x9E3779B185EBCA87)) % mod

  def _guard_append(T, items, cap):
    if not items:
      return T
    room = cap - len(T)
    if room <= 0:
      return T
    if len(items) > room:
      items = items[:room]
    T.extend(items)
    return T

  def _append_connectors(S, starts, delta):
    # Classic four connectors help enforce FF pressure without blowing up omega.
    s0, s1, s2, s3 = starts
    S.append(((s0 - 1) * delta, (s1 - 1) * delta))  # left cap
    S.append(((s2 + 2) * delta, (s3 + 2) * delta))  # right cap
    S.append(((s0 + 2) * delta, (s2 - 1) * delta))  # cross 1
    S.append(((s1 + 2) * delta, (s3 - 1) * delta))  # cross 2

  def _apply_round(current_T, starts, ridx, always_interleave=True):
    lo, hi, delta = _span_delta(current_T)

    # Build four translated blocks
    blocks = []
    for s in starts:
      base = s * delta - lo
      block = [(l + base, r + base) for (l, r) in current_T]
      blocks.append(block)

    # Interleaving every round; block order seed-derived for diversity
    S = []
    if always_interleave:
      order = list(range(4))
      if (_det_choice(BASE_SEED, ridx, 2) == 1):
        order.reverse()
      maxlen = max(len(b) for b in blocks)
      for i in range(maxlen):
        for idx in order:
          blk = blocks[idx]
          if i < len(blk):
            S.append(blk[i])
    else:
      # Fallback: sequential blocks, seed-derived reverse
      block_seq = list(blocks)
      if (_det_choice(BASE_SEED, ridx, 2) == 1):
        block_seq.reverse()
      for blk in block_seq:
        S.extend(blk)

    # Append classic connectors
    _append_connectors(S, starts, delta)
    return S

  # Small density pins inserted occasionally at seed-derived positions
  def _spine_pins(current_T, per_round_quota=0, ridx=0):
    if per_round_quota <= 0 or not current_T:
      return []
    lo, hi, delta = _span_delta(current_T)
    G = max(1, hi - lo)
    eps = max(1, G // 1024)
    # choose 2 fractional anchors per round deterministically
    anchors = []
    a1 = (12 + _det_choice(BASE_SEED, ridx * 3 + 1, 8)) / 100.0  # ~[0.12..0.19]
    a2 = (72 + _det_choice(BASE_SEED, ridx * 3 + 2, 8)) / 100.0  # ~[0.72..0.79]
    anchors.extend([a1, a2])
    pins = []
    for j in range(per_round_quota):
      f = anchors[j % len(anchors)]
      L = lo + int(round(f * G)) + (j % 3)
      R = L + eps
      if R > L:
        pins.append((L, R))
    return pins

  # -------------------- Stage 0: seed --------------------
  if seed_count <= 1:
    T = [(0, 1)]
  else:
    step = 3
    T = [(i * step, i * step + 1) for i in range(seed_count)]

  # -------------------- Stage 1: six-template, seed-driven interleaved spine --------------------
  # Guard the final spine round so we keep SPINE_RESERVE slots for micro phases.
  pin_allowance = SPINE_PIN_BUDGET
  for ridx in range(SPINE_ROUNDS_MAX):
    if _predict_next_size(len(T)) > CAP - SPINE_RESERVE:
      break
    # Seed-driven template selection
    t_idx = _det_choice(BASE_SEED, ridx, len(template_bank))
    starts = template_bank[t_idx]
    T = _apply_round(T, starts, ridx, always_interleave=SPINE_ALWAYS_INTERLEAVE)
    # Inject a few short pins after this round (tiny count, CAP-guarded)
    if pin_allowance > 0:
      per_round = 2 if pin_allowance >= 2 else 1
      pins = _spine_pins(T, per_round_quota=per_round, ridx=ridx)
      T = _guard_append(T, pins, CAP)
      pin_allowance -= per_round

  if len(T) >= CAP:
    return T[:CAP]

  # -------------------- Stage 1.5: deterministic long-range connectors (CAP-guarded) --------------------
  if LR_CONNECTORS_ENABLED and LR_CONNECTOR_BUDGET > 0:
    lo, hi, delta = _span_delta(T)
    lr = []
    # two or three conservative cross-scale connectors
    lr.append((lo + (2 + 4) * delta, lo + (14 + 4) * delta))
    lr.append((lo + (6 + 5) * delta, lo + (10 + 5) * delta))
    lr.append((lo + (3 + 6) * delta, lo + (15 + 6) * delta))
    lr = [(a, b) for (a, b) in lr if b > a]
    lr = lr[:LR_CONNECTOR_BUDGET]
    T = _guard_append(T, lr, CAP)

  if len(T) >= CAP:
    return T[:CAP]

  # -------------------- Stage 2: dual micro-phase with explicit budgets --------------------
  def _thin_seed(current_T, max_seed):
    n = len(current_T)
    if n == 0 or max_seed <= 0:
      return []
    step = max(1, n // max_seed)
    return current_T[::step][:max_seed]

  def _build_micro(current_T, budget, family_id=0, iter_id=0):
    if not current_T or budget <= 0:
      return []
    glo = min(l for l, r in current_T)
    ghi = max(r for l, r in current_T)
    G = max(1, ghi - glo)

    # Thin evenly spaced seed, slightly different sizing per family
    seed_sz = max(8, min(36 if family_id == 0 else 40, len(current_T) // (280 if family_id == 0 else 240)))
    U = _thin_seed(current_T, seed_sz)
    if not U:
      return []

    ulo = min(l for l, r in U)

    # Window families with tiny deterministic shifts
    if family_id == 0:
      base_windows = [(0.12, 0.22), (0.36, 0.46)]
    else:
      base_windows = [(0.62, 0.72), (0.84, 0.94)]

    def _shift(frac, k):
      # tiny de-aliasing shift per (family_id, iter_id, window index)
      raw = _det_choice(BASE_SEED ^ (family_id << 8) ^ (iter_id << 4), k, 11)
      sign = -1 if (raw % 2 == 0) else 1
      mag = (raw % 5) * 0.001  # up to 0.004
      return max(0.03, min(0.95, frac + sign * mag))

    windows = []
    for k, (fa, fb) in enumerate(base_windows):
      a = _shift(fa, k)
      b = _shift(fb, k + 1)
      if b <= a + 0.01:
        b = a + 0.02
      b = min(0.97, b)
      windows.append((a, b))

    # Place translated micro-blocks aligned to windows
    blocks = []
    for k, (fa, fb) in enumerate(windows):
      win_lo = glo + int(round(fa * G))
      base = win_lo - ulo
      block = [(l + base, r + base) for (l, r) in U]
      # alternate internal reversal
      if ((k + iter_id + family_id) % 2) == 1:
        block = list(reversed(block))
      blocks.append(block)

    # Interleave micro blocks
    micro = []
    maxlen = max(len(b) for b in blocks) if blocks else 0
    order = list(range(len(blocks)))
    if (family_id + iter_id) % 2 == 1:
      order.reverse()
    for i in range(maxlen):
      for idx in order:
        blk = blocks[idx]
        if i < len(blk):
          micro.append(blk[i])

    # Fractional-span connectors (family-dependent)
    connectors = []
    if family_id == 0:
      connectors.append((glo + int(round(0.10 * G)), glo + int(round(0.30 * G))))
      connectors.append((glo + int(round(0.32 * G)), glo + int(round(0.52 * G))))
    else:
      connectors.append((glo + int(round(0.60 * G)), glo + int(round(0.82 * G))))
      connectors.append((glo + int(round(0.78 * G)), glo + int(round(0.96 * G))))
    micro.extend([(a, b) for a, b in connectors if b > a])

    # Trim to budget
    if len(micro) > budget:
      micro = micro[:budget]
    return micro

  # Compute budgets explicitly
  room_after_spine = max(0, CAP - len(T))
  if room_after_spine > 0:
    budget_A = int(room_after_spine * PHASE2_BUDGET_FRACTION_A)
    budget_B = room_after_spine - budget_A
  else:
    budget_A = budget_B = 0

  # Micro-phase A
  if budget_A > 0:
    microA = _build_micro(T, budget_A, family_id=0, iter_id=0)
    T = _guard_append(T, microA, CAP)

  # Micro-phase B
  if budget_B > 0:
    microB = _build_micro(T, budget_B, family_id=1, iter_id=1)
    T = _guard_append(T, microB, CAP)

  # -------------------- Stage 3: tail connectors (small, CAP-guarded) --------------------
  if LR_CONNECTORS_ENABLED and len(T) < CAP:
    lo, hi, delta = _span_delta(T)
    G = max(1, hi - lo)
    # three conservative tail connectors
    tails = [
      (lo + int(round(0.08 * G)), lo + int(round(0.58 * G))),
      (lo + int(round(0.42 * G)), lo + int(round(0.72 * G))),
      (lo + int(round(0.76 * G)), lo + int(round(0.94 * G))),
    ]
    tails = [(a, b) for (a, b) in tails if b > a]
    T = _guard_append(T, tails[:max(0, CAP - len(T))], CAP)

  # Final CAP trim
  if len(T) > CAP:
    T = T[:CAP]
  return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()