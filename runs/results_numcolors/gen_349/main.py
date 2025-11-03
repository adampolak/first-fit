# EVOLVE-BLOCK-START

def construct_intervals(seed_count=1):
  """
  Deterministic 6-template KT-style backbone with two-pass micro phases and
  deterministic long-range connectors. Returns a list of open intervals (l, r)
  in FirstFit presentation order.

  Inputs:
    seed_count (int): preserved for interface compatibility (unused).

  Output:
    List[(int, int)]: open intervals in FF order.
  """

  # Global capacity to keep total intervals < 10000
  CAP = 9800

  # Deterministic base seed for template and connector choices
  BASE_SEED = 911

  # Backbone depth (KT rounds)
  SPINE_ROUNDS = 6

  # Interleaving policy for spine
  INTERLEAVE_ON_EVEN = True

  # Six deterministic start templates (diversified offsets)
  TEMPLATE_BANK = [
    (2, 6, 10, 14),  # T1 classic KT
    (1, 5, 9, 13),   # T2 left-shifted
    (3, 7, 11, 15),  # T3 right-shifted
    (4, 8, 12, 16),  # T4 stretched-right
    (2, 5, 11, 14),  # T5 skew-tight
    (1, 7, 11, 15),  # T6 wide skew
  ]

  # Seed with a single unit interval to avoid early omega inflation
  T = [(0, 1)]

  # ============== Helpers ==============

  def _span_delta(current_T):
    lo = min(l for l, r in current_T)
    hi = max(r for l, r in current_T)
    delta = hi - lo
    if delta <= 0:
      delta = 1
    return lo, hi, delta

  def _normalize(current_T):
    """Translate so that min l is 0. Preserves all intersections."""
    if not current_T:
      return current_T
    lo = min(l for l, _ in current_T)
    if lo == 0:
      return current_T
    return [(l - lo, r - lo) for (l, r) in current_T]

  def _append_connectors(S, starts, delta):
    # Classic four connectors; strong FF pressure without obvious omega blow-up.
    s0, s1, s2, s3 = starts
    S.append(((s0 - 1) * delta, (s1 - 1) * delta))  # left cap
    S.append(((s2 + 2) * delta, (s3 + 2) * delta))  # right cap
    S.append(((s0 + 2) * delta, (s2 - 1) * delta))  # cross 1
    S.append(((s1 + 2) * delta, (s3 - 1) * delta))  # cross 2

  def _apply_spine_round(current_T, starts, ridx, interleave_even=True):
    lo, hi, delta = _span_delta(current_T)

    # Build four translated blocks
    blocks = []
    for s in starts:
      base = s * delta - lo
      block = [(l + base, r + base) for (l, r) in current_T]
      blocks.append(block)

    # Controlled interleaving: on even rounds; reverse order on odd
    do_interleave = (interleave_even and (ridx % 2 == 0))
    reverse_order = (ridx % 2 == 1)

    S = []
    if do_interleave:
      order = list(range(4))
      if reverse_order:
        order.reverse()
      maxlen = max(len(b) for b in blocks)
      for i in range(maxlen):
        for idx in order:
          blk = blocks[idx]
          if i < len(blk):
            S.append(blk[i])
    else:
      if reverse_order:
        blocks = list(reversed(blocks))
      for blk in blocks:
        S.extend(blk)

    # Append classic connectors
    _append_connectors(S, starts, delta)
    return S

  # ============== Stage 1: Six-round KT backbone with six-template rotation ==============

  for ridx in range(SPINE_ROUNDS):
    # Predict next size: |T_next| ~ 4*|T| + 4; guard CAP
    if 4 * len(T) + 4 > CAP:
      break
    # Deterministic template selection
    tidx = (BASE_SEED + 7 * ridx + ridx * ridx) % len(TEMPLATE_BANK)
    starts = TEMPLATE_BANK[tidx]
    T = _apply_spine_round(T, starts, ridx, interleave_even=INTERLEAVE_ON_EVEN)
    if len(T) >= CAP:
      T = T[:CAP]
      return T

  T = _normalize(T)

  # Early guard
  if len(T) >= CAP - 8:
    return T

  # ============== Stage 1.5: Deterministic long-range connectors ==============
  # 4–6 connectors across fractional span; conservative to control omega.
  lo, hi, delta = _span_delta(T)
  G = max(1, hi - lo)

  def _frac_iv(a, b):
    L = lo + max(1, int(round(a * G)))
    R = lo + max(1, int(round(b * G)))
    if R <= L:
      R = L + 1
    return (L, R)

  # Connector schedule derived from BASE_SEED to break symmetry consistently
  LR_CONNECTORS = [
    _frac_iv(0.06, 0.33),
    _frac_iv(0.18, 0.52),
    _frac_iv(0.34, 0.79),
    _frac_iv(0.57, 0.90),
    _frac_iv(0.12, 0.88),
    _frac_iv(0.41, 0.69),
  ]

  # Gate connector count by remaining capacity
  room = CAP - len(T)
  if room > 0:
    add = LR_CONNECTORS[:min(len(LR_CONNECTORS), room)]
    T.extend(add)

  if len(T) >= CAP - 16:
    return _normalize(T)

  # ============== Stage 2: Two-pass micro phases with adaptive budgets ==============

  # Budget split: two-phase micro injection; deterministic and CAP-aware
  remaining = max(0, CAP - len(T))
  budgetA_frac = 0.62
  budgetA = int(remaining * budgetA_frac)
  budgetB = max(0, remaining - budgetA)

  def _build_micro_round(current_T, budget, iter_id=0, alt=False):
    # Defensive guards
    if not current_T or budget <= 8:
      return []

    glo = min(l for l, r in current_T)
    ghi = max(r for l, r in current_T)
    G = max(1, ghi - glo)

    # Thin, evenly spaced seed; tie to budget to stay lean
    # Increase density slightly on alt phase to strengthen coupling
    base_seed_cap = 36 if not alt else 44
    seed_sz = max(8, min(base_seed_cap, max(8, budget // 32)))
    stride = max(1, len(current_T) // seed_sz)
    U = [current_T[i] for i in range(0, len(current_T), stride)][:seed_sz]
    if not U:
      return []

    ulo = min(l for l, r in U)

    # Two distinct window families A/B; B has slight inward shifts and different coverage
    if not alt:
      # Phase A windows (five micro windows, spread inside [0.06, 0.92])
      window_fracs = [
        (0.10, 0.18),
        (0.28, 0.36),
        (0.46, 0.54),
        (0.64, 0.72),
        (0.82, 0.90),
      ]
      interleave_reverse = (iter_id % 2 == 1)
    else:
      # Phase B windows (shifted, asymmetric to catch missed couplings)
      window_fracs = [
        (0.06, 0.14),
        (0.24, 0.33),
        (0.40, 0.49),
        (0.58, 0.66),
        (0.76, 0.85),
      ]
      interleave_reverse = (iter_id % 2 == 0)

    # Build translated micro-blocks aligned to these windows
    blocks = []
    for (fa, fb) in window_fracs:
      win_lo = glo + int(round(fa * G))
      base = win_lo - ulo
      block = [(l + base, r + base) for (l, r) in U]
      # Deterministic internal reversal by window to break symmetry
      tag = (int(round(fa * 100)) // 5 + iter_id + (1 if alt else 0)) % 2
      if tag == 1:
        block = list(reversed(block))
      blocks.append(block)

    # Interleave micro-blocks
    micro = []
    maxlen = max(len(b) for b in blocks)
    order = list(range(len(blocks)))
    if interleave_reverse:
      order.reverse()
    for i in range(maxlen):
      for idx in order:
        blk = blocks[idx]
        if i < len(blk):
          micro.append(blk[i])

    # Deterministic fractional-span connectors, conservative coverage
    micro_connectors = []
    if not alt:
      # A: longer pair + two mid crosses
      micro_connectors.extend([
        (glo + int(round(0.08 * G)), glo + int(round(0.30 * G))),
        (glo + int(round(0.62 * G)), glo + int(round(0.92 * G))),
        (glo + int(round(0.26 * G)), glo + int(round(0.56 * G))),
        (glo + int(round(0.44 * G)), glo + int(round(0.78 * G))),
      ])
    else:
      # B: slightly shorter pair + one long-range coupler
      micro_connectors.extend([
        (glo + int(round(0.12 * G)), glo + int(round(0.44 * G))),
        (glo + int(round(0.55 * G)), glo + int(round(0.86 * G))),
        (glo + int(round(0.18 * G)), glo + int(round(0.84 * G))),
      ])

    for a, b in micro_connectors:
      if b > a:
        micro.append((a, b))

    # Enforce budget strictly
    if len(micro) > budget:
      micro = micro[:budget]
    return micro

  # Phase A
  if budgetA > 0:
    microA = _build_micro_round(T, budgetA, iter_id=0, alt=False)
    if microA:
      room = CAP - len(T)
      if len(microA) > room:
        microA = microA[:room]
      T.extend(microA)

  # Phase B (recompute remaining)
  budgetB = max(0, CAP - len(T))
  if budgetB > 0:
    microB = _build_micro_round(T, budgetB, iter_id=1, alt=True)
    if microB:
      room = CAP - len(T)
      if len(microB) > room:
        microB = microB[:room]
      T.extend(microB)

  # Final normalization and capacity trim
  if len(T) > CAP:
    T = T[:CAP]
  T = _normalize(T)
  return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()