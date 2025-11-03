# EVOLVE-BLOCK-START

def construct_intervals(seed_count=1):
  """
  Enhanced deterministic spine + adaptive micro-phases generator.

  Interface preserved: construct_intervals(seed_count=1) -> list[(l, r)]
  """

  # Capacity guard (strict)
  CAP = 9800

  # Expanded deterministic six-template bank (diverse relative offsets)
  TEMPLATE_BANK = [
    (2, 6, 10, 14),  # classic KT
    (1, 5, 9, 13),   # left shift
    (3, 7, 11, 15),  # right shift
    (4, 8, 12, 16),  # stretched-right
    (2, 5, 8, 11),   # compressed stepping
    (3, 6, 9, 12),   # alternate compressed
  ]

  # Deterministic interleave-round pattern (apply interleaving on these indices)
  INTERLEAVE_PATTERN = {1, 2, 4}  # interleave on rounds 1,2,4 (0-based indexing)

  # Seed: single unit interval to avoid early omega inflation.
  T = [(0, 1)] if seed_count <= 1 else [(0, 1)]

  # ---------- Helpers ----------
  def _span_delta(current_T):
    # Single pass min/max for speed and robustness
    if not current_T:
      return 0, 0, 1
    lo = current_T[0][0]
    hi = current_T[0][1]
    for l, r in current_T:
      if l < lo: lo = l
      if r > hi: hi = r
    delta = hi - lo
    if delta <= 0:
      delta = 1
    return lo, hi, delta

  def _rebase_to_zero(current_T):
    if not current_T:
      return current_T
    lo, _, _ = _span_delta(current_T)
    if lo == 0:
      return current_T
    return [(l - lo, r - lo) for (l, r) in current_T]

  def _append_connectors(S, starts, delta):
    # Classic four connectors kept plus modest long-range connector candidate
    s0, s1, s2, s3 = starts
    S.append(((s0 - 1) * delta, (s1 - 1) * delta))  # left cap
    S.append(((s2 + 2) * delta, (s3 + 2) * delta))  # right cap
    S.append(((s0 + 2) * delta, (s2 - 1) * delta))  # cross 1
    S.append(((s1 + 2) * delta, (s3 - 1) * delta))  # cross 2

  def _cap_at(lo, span, a_frac, b_frac):
    L = lo + max(1, int(round(a_frac * span)))
    R = lo + max(1, int(round(b_frac * span)))
    if R <= L:
      R = L + 1
    return (L, R)

  def _apply_round(current_T, starts, do_interleave=False, reverse_order=False, density_mult=1.0):
    lo, hi, delta = _span_delta(current_T)

    # Build translated blocks
    base_vals = [s * delta - lo for s in starts]
    blocks = [[(l + base, r + base) for (l, r) in current_T] for base in base_vals]

    # Assemble sequence S either interleaving or sequential
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

    # Append connectors
    _append_connectors(S, starts, delta)

    # CAP-aware density multiplier:
    # If density_mult > 1.0, duplicate a tiny deterministic sample of S (spread by small offsets)
    # but do not exceed CAP and do not grow clique significantly (use tiny offsets inside delta).
    if density_mult > 1.001 and len(S) < CAP:
      # Choose deterministic sample size proportional to S and density multiplier
      extra_budget = int(max(1, min(CAP - len(S), int((density_mult - 1.0) * len(S) / 2))))
      # Select every k-th element deterministically
      if extra_budget > 0 and len(S) > 0:
        step = max(1, len(S) // extra_budget)
        # small offset fraction to avoid making intervals identical and avoid large cliques
        offset = max(1, delta // 100)  # at least 1
        # Insert duplicates near the middle to pressure FF where many colors are active
        mid = len(S) // 2
        inserts = []
        for i in range(0, len(S), step):
          if len(inserts) >= extra_budget:
            break
          l, r = S[i]
          # Shift by a deterministic sign pattern to spread but keep intersections similar
          sign = 1 if ((i // step) % 2 == 0) else -1
          nl = l + sign * offset
          nr = r + sign * offset
          if nr <= nl:
            nr = nl + 1
          inserts.append((nl, nr))
        # Insert near mid as alternating insertions to maximize overlap with active colors
        pos = max(0, mid - len(inserts) // 2)
        for j, iv in enumerate(inserts):
          if len(S) >= CAP:
            break
          S.insert(min(len(S), pos + 2 * j), iv)

    return S

  # Predictive cap-aware number of rounds (4n + 4 growth)
  def _max_rounds_within_cap(initial_size, max_rounds):
    sz = initial_size
    done = 0
    for _ in range(max(0, int(max_rounds))):
      nxt = 4 * sz + 4
      if nxt > CAP:
        break
      sz = nxt
      done += 1
    return done

  # ---------- Stage 1: KT-style spine with six-template rotation ----------
  target_rounds = _max_rounds_within_cap(len(T), 6)
  # Use a slight adaptive density multiplier schedule per round (bounded)
  density_schedule = [1.00, 1.06, 1.04, 1.08, 1.03, 1.05]  # deterministic multipliers

  for ridx in range(target_rounds):
    starts = TEMPLATE_BANK[ridx % len(TEMPLATE_BANK)]
    # Deterministic interleaving rule from pattern
    do_inter = (ridx in INTERLEAVE_PATTERN)
    rev = (ridx % 2 == 1)
    # Use density multiplier from schedule (cycle if fewer entries)
    density_mult = density_schedule[ridx % len(density_schedule)]
    # Apply round
    T = _apply_round(T, starts, do_interleave=do_inter, reverse_order=rev, density_mult=density_mult)
    # Rebase to keep numbers stable
    if len(T) >= CAP:
      T = T[:CAP]
      return T
    T = _rebase_to_zero(T)

  # Early return guard
  if len(T) >= CAP - 12:
    return T

  # ---------- Deterministic long-range connectors appended after spine ----------
  lo, hi, span = _span_delta(T)
  # Add a few long connectors that span large fractions of the current span.
  long_connectors = [
    _cap_at(lo, span, 0.02, 0.95),
    _cap_at(lo, span, 0.10, 0.80),
    _cap_at(lo, span, 0.30, 0.90),
  ]
  for c in long_connectors:
    if len(T) >= CAP:
      break
    # Append near the tail to maximize overlap with active colors
    T.insert(len(T) - 3 if len(T) >= 3 else len(T), c)
  T = _rebase_to_zero(T)

  if len(T) >= CAP - 12:
    return T

  # ---------- Adaptive two-phase micro budgeting ----------
  remaining = CAP - len(T)
  # Phase A: allocate a fraction (e.g., 42%) of remaining to primary micro windows
  phaseA_budget = max(0, int(round(remaining * 0.42)))
  # Reserve small safety headroom
  if phaseA_budget > remaining - 8:
    phaseA_budget = max(0, remaining - 8)

  # Phase B: leftover budget (but keep at least 4 slots for final caps)
  remaining_after_A = CAP - len(T) - phaseA_budget
  phaseB_budget = max(0, min(remaining_after_A - 4, CAP - len(T) - phaseA_budget))
  # Safety cap
  phaseA_budget = min(phaseA_budget, 8000)
  phaseB_budget = min(phaseB_budget, 8000)

  def _thin_seed(current_T, max_seed):
    n = len(current_T)
    if n == 0 or max_seed <= 0:
      return []
    step = max(1, n // max_seed)
    return current_T[::step][:max_seed]

  def _build_micro(current_T, budget, iter_id=0, alt=False):
    if not current_T or budget <= 4:
      return []

    glo, ghi, G = _span_delta(current_T)
    # Seed: adaptively sized but capped
    seed_sz = max(8, min(48, len(current_T) // 200))
    U = _thin_seed(current_T, seed_sz)
    if not U:
      return []

    ulo = min(l for l, _ in U)

    # Two families of windows; alt uses offset windows to catch missed pairings
    if not alt:
      shift = ((iter_id % 4) * 0.015)
      window_fracs = [
        (0.10 + shift, 0.20 + shift),
        (0.32 + shift, 0.42 + shift),
        (0.54 + shift, 0.64 + shift),
        (0.76 + shift, 0.86 + shift),
      ]
    else:
      # alternate windows slightly overlapping different spine parts
      window_fracs = [
        (0.05, 0.15),
        (0.28, 0.38),
        (0.50, 0.62),
        (0.70, 0.84),
      ]

    blocks = []
    for idx, (fa, fb) in enumerate(window_fracs):
      win_lo = glo + int(round(max(0.0, min(0.99, fa)) * G))
      base = win_lo - ulo
      block = [(l + base, r + base) for (l, r) in U]
      # Deterministic internal reversal to disrupt symmetry
      if (idx + iter_id) % 2 == 1:
        block = list(reversed(block))
      blocks.append(block)

    # Interleave blocks to maximize cross-color hits
    micro = []
    maxlen = max(len(b) for b in blocks)
    order = list(range(len(blocks)))
    if iter_id % 2 == 1:
      order.reverse()
    for i in range(maxlen):
      for j in order:
        blk = blocks[j]
        if i < len(blk):
          micro.append(blk[i])

    # Micro connectors at fractional scale (careful to not inflate omega)
    micro_connectors = [
      (glo + int(round(0.06 * G)), glo + int(round(0.28 * G))),
      (glo + int(round(0.34 * G)), glo + int(round(0.58 * G))),
      (glo + int(round(0.62 * G)), glo + int(round(0.90 * G))),
    ]
    if alt:
      # add one more cross in alternate
      micro_connectors.append((glo + int(round(0.18 * G)), glo + int(round(0.82 * G))))

    for a, b in micro_connectors:
      if b > a:
        micro.append((a, b))

    # Trim micro to budget
    if len(micro) > budget:
      micro = micro[:budget]
    return micro

  # Execute Phase A
  if phaseA_budget >= 8:
    microA = _build_micro(T, phaseA_budget, iter_id=0, alt=False)
    if microA:
      # Insert microA evenly across the tail region to maximize simultaneous active colors
      insert_stride = max(1, len(microA) // 6)
      pos = len(T) - 2
      for i, iv in enumerate(microA):
        if len(T) >= CAP:
          break
        # alternate insertion positions to overlap active color windows
        T.insert(max(0, pos - (i % 3)), iv)
      T = _rebase_to_zero(T)

  # Execute Phase B (alternate windows) with leftover budget
  remaining = CAP - len(T)
  if remaining <= 8:
    # small final near-tail caps if anything remains
    if remaining > 0:
      lo, hi, sp = _span_delta(T)
      final_caps = [
        _cap_at(lo, sp, 0.12, 0.55),
        _cap_at(lo, sp, 0.40, 0.85)
      ]
      for cap in final_caps[:remaining]:
        T.append(cap)
    T = _rebase_to_zero(T)
    if len(T) > CAP:
      T = T[:CAP]
    return T

  phaseB_budget = remaining - 4 if remaining - 4 >= 8 else remaining - 1
  if phaseB_budget >= 6:
    microB = _build_micro(T, phaseB_budget, iter_id=1, alt=True)
    if microB:
      # Append microB primarily near the end, mixing with near-tail caps
      for i, iv in enumerate(microB):
        if len(T) >= CAP:
          break
        # place microB items near end with slight stagger
        pos = len(T) - 1 - (i % 5)
        if pos < 0:
          T.append(iv)
        else:
          T.insert(pos, iv)
      T = _rebase_to_zero(T)

  # Final small near-tail caps to use leftover budget and tie colors across windows
  remaining = CAP - len(T)
  if remaining > 0:
    lo, hi, sp = _span_delta(T)
    tail_caps = [
      _cap_at(lo, sp, 0.08, 0.60),
      _cap_at(lo, sp, 0.25, 0.75),
      _cap_at(lo, sp, 0.75, 0.92),
    ]
    for cap in tail_caps[:remaining]:
      T.append(cap)
    T = _rebase_to_zero(T)

  # Final trim safety
  if len(T) > CAP:
    T = T[:CAP]

  return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()