# EVOLVE-BLOCK-START

def construct_intervals(enable_alt_microphase=True):
  """
  Construct a sequence of intervals of real line,
  in the order presented to FirstFit, to maximize FF colors / omega,
  while keeping omega around <= 10.

  Novel elements:
  - Six-template KT spine with deterministic interleaving schedule.
  - Two-pass CAP-aware micro phases anchored to spine span (deterministic windows).
  - Long-range connector augmentation pass after spine.
  - Deterministic central density boost after each spine round.
  - Tower + caps gadget to couple many colors across layers at low omega cost.

  Args:
    enable_alt_microphase (bool): enable the second window family micro-phase.

  Returns:
    list[(int,int)]: open intervals (l, r) in the arrival order for FirstFit.
  """

  # Hard capacity guard
  CAP = 9800

  # Six-template rotation (deterministic)
  template_bank = [
    (2, 6, 10, 14),  # T1 classic
    (1, 5, 9, 13),   # T2 left-shift
    (3, 7, 11, 15),  # T3 right-shift
    (4, 8, 12, 16),  # T4 stretched-right
    (5, 9, 13, 17),  # T5 extra right-shift
    (6, 10, 14, 18), # T6 extra classic-shift
  ]

  # Seed: single unit interval keeps early omega tiny
  T = [(0, 1)]

  # Helpers
  def _span_of(seq):
    lo = min(l for l, r in seq)
    hi = max(r for l, r in seq)
    delta = hi - lo
    if delta <= 0:
      delta = 1
    return lo, hi, int(delta)

  def _apply_spine_round(current_T, starts, do_interleave=False, reverse_order=False):
    lo, hi, delta = _span_of(current_T)

    # Build blocks (four translated copies)
    blocks = []
    for s in starts:
      base = s * delta - lo
      blk = [(l + base, r + base) for (l, r) in current_T]
      blocks.append(blk)

    # Interleaving policy (deterministic, parity-aware)
    S = []
    if do_interleave:
      # Order by start value; reverse when requested
      order = sorted(range(4), key=lambda i: starts[i], reverse=reverse_order)
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

    # Classic KT connectors
    s0, s1, s2, s3 = starts
    S.append(((s0 - 1) * delta + lo, (s1 - 1) * delta + lo))  # left cap
    S.append(((s2 + 2) * delta + lo, (s3 + 2) * delta + lo))  # right cap
    S.append(((s0 + 2) * delta + lo, (s2 - 1) * delta + lo))  # cross 1
    S.append(((s1 + 2) * delta + lo, (s3 - 1) * delta + lo))  # cross 2
    return S

  def _central_density_boost(seq, budget_per_round=6):
    # Add very short pins near the central span to raise FF pressure slightly.
    if not seq or budget_per_round <= 0:
      return []
    lo, hi, delta = _span_of(seq)
    G = max(1, hi - lo)
    eps = max(1, G // 4096)
    base = lo + G // 2
    stride = max(3, 2 * eps + 3)
    pins = []
    start = base - (budget_per_round // 2) * stride
    for i in range(budget_per_round):
      L = start + i * stride
      R = L + eps
      if R > L:
        pins.append((L, R))
    return pins

  # Predictive guard for rounds: size -> 4*size + 4 per round
  def _max_rounds(sz, max_rounds):
    done = 0
    for _ in range(max_rounds):
      nxt = 4 * sz + 4
      if nxt > CAP:
        break
      sz = nxt
      done += 1
    return done

  # Stage 1: Six-round spine with deterministic interleaving and parity reversal
  rounds = _max_rounds(len(T), 6)
  for ridx in range(rounds):
    starts = template_bank[ridx % len(template_bank)]
    do_inter = (ridx % 2 == 0)     # interleave on even rounds
    rev = (ridx % 2 == 1)          # reverse block order on odd rounds
    T = _apply_spine_round(T, starts, do_interleave=do_inter, reverse_order=rev)
    # Small density boost after each round (CAP-guarded)
    if len(T) < CAP - 8:
      pins = _central_density_boost(T, budget_per_round=6 if ridx < 4 else 4)
      room = CAP - len(T)
      if pins and room > 0:
        if len(pins) > room:
          pins = pins[:room]
        T.extend(pins)
    if len(T) >= CAP:
      T = T[:CAP]
      return T

  if len(T) >= CAP - 64:
    return T[:CAP]

  # Long-range connector augmentation pass (6 connectors)
  def _augment_connectors(seq, max_cons=6):
    if not seq or max_cons <= 0:
      return []
    lo, hi, delta = _span_of(seq)
    G = max(1, hi - lo)
    cons = []
    # Fractional long caps and crosses; lengths near broad fractions of G.
    fracs = [
      (0.08, 0.32), (0.70, 0.92),  # two caps
      (0.18, 0.84), (0.26, 0.74),  # two crosses
      (0.35, 0.65), (0.42, 0.58),  # mid connectors
    ]
    for (a, b) in fracs[:max_cons]:
      A = lo + int(round(a * G))
      B = lo + int(round(b * G))
      if B > A:
        cons.append((A, B))
    return cons

  aug = _augment_connectors(T, max_cons=6)
  if aug:
    room = CAP - len(T)
    if len(aug) > room:
      aug = aug[:room]
    T.extend(aug)
  if len(T) >= CAP - 64:
    return T[:CAP]

  # Two-pass CAP-aware micro phases anchored to the spine span
  # Phase 1 uses primary windows; Phase 2 (optional) uses alternate windows.
  def _thin_seed(current_T, max_seed):
    n = len(current_T)
    if n <= 0 or max_seed <= 0:
      return []
    step = max(1, n // max_seed)
    return current_T[::step][:max_seed]

  def _build_window_micro(seq, budget, windows, reverse=False, add_connectors=True):
    if not seq or budget <= 8:
      return []
    glo, ghi, G = _span_of(seq)
    # Keep seed small to limit omega growth
    seed_sz = max(10, min(40, len(seq) // 250))
    U = _thin_seed(seq, seed_sz)
    if not U:
      return []
    ulo = min(l for l, r in U)
    blocks = []
    # Create blocks in each window by translating the seed
    for (fa, fb) in windows:
      win_lo = glo + int(round(fa * G))
      base = win_lo - ulo
      blk = [(l + base, r + base) for (l, r) in U]
      if reverse:
        blk = list(reversed(blk))
      blocks.append(blk)

    # Interleave blocks deterministically
    micro = []
    maxlen = max(len(b) for b in blocks)
    order = list(range(len(blocks)))
    if reverse:
      order.reverse()
    for i in range(maxlen):
      for idx in order:
        blk = blocks[idx]
        if i < len(blk):
          micro.append(blk[i])

    # Add fractional connectors across the windows (optional)
    if add_connectors:
      cons = [
        (glo + int(round(0.08 * G)), glo + int(round(0.30 * G))),
        (glo + int(round(0.60 * G)), glo + int(round(0.92 * G))),
        (glo + int(round(0.26 * G)), glo + int(round(0.56 * G))),
        (glo + int(round(0.44 * G)), glo + int(round(0.78 * G))),
      ]
      for (a, b) in cons:
        if b > a:
          micro.append((a, b))

    if len(micro) > budget:
      micro = micro[:budget]
    return micro

  # Phase budgets: deterministic split of remaining CAP
  remaining = CAP - len(T)
  if remaining > 24:
    phase1_budget = max(0, int(0.55 * remaining))
    phase2_budget = max(0, remaining - phase1_budget)

    primary_windows = [(0.12, 0.22), (0.35, 0.45), (0.58, 0.68), (0.80, 0.90)]
    micro1 = _build_window_micro(T, phase1_budget, primary_windows, reverse=False, add_connectors=True)
    if micro1:
      room = CAP - len(T)
      if len(micro1) > room:
        micro1 = micro1[:room]
      T.extend(micro1)

    # Optional Phase 2 with alternate windows (smaller budget improves stability)
    if enable_alt_microphase:
      remaining = CAP - len(T)
      if remaining > 12 and phase2_budget > 0:
        alt_windows = [(0.08, 0.18), (0.28, 0.38), (0.62, 0.72), (0.82, 0.92)]
        phase2_budget = min(phase2_budget, remaining)
        micro2 = _build_window_micro(T, phase2_budget, alt_windows, reverse=True, add_connectors=True)
        if micro2:
          room = CAP - len(T)
          if len(micro2) > room:
            micro2 = micro2[:room]
          T.extend(micro2)

  if len(T) >= CAP - 64:
    return T[:CAP]

  # Tower + caps gadget: three localized stacks with staggered short intervals.
  def _build_towers(seq, budget, towers=3, layers=10):
    if not seq or budget <= towers:
      return []
    glo, ghi, G = _span_of(seq)
    eps = max(1, G // 4096)
    # Tower centers in separated windows to avoid clique spikes
    anchors = [0.18, 0.46, 0.74][:towers]
    stride = max(3, 2 * eps + 3)
    out = []
    for ti, frac in enumerate(anchors):
      base = glo + int(round(frac * G))
      start = base - (layers // 2) * stride + (ti % 2)
      for j in range(layers):
        L = start + j * stride
        R = L + eps
        if R > L:
          out.append((L, R))
    # Caps that touch one piece from each tower
    caps = []
    if towers >= 3:
      caps.append((glo + int(round(0.17 * G)), glo + int(round(0.75 * G))))
      caps.append((glo + int(round(0.20 * G)), glo + int(round(0.52 * G))))
      caps.append((glo + int(round(0.48 * G)), glo + int(round(0.82 * G))))
    out.extend(caps)
    return out[:budget]

  room = CAP - len(T)
  if room > 8:
    # Keep towers modest to avoid omega blow-up; gate by remaining room.
    layers = min(18, max(8, room // 300))
    towers = 3
    tower_pack = _build_towers(T, budget=room, towers=towers, layers=layers)
    if tower_pack:
      if len(tower_pack) > room:
        tower_pack = tower_pack[:room]
      T.extend(tower_pack)

  # Final safety: trim to CAP
  if len(T) > CAP:
    T = T[:CAP]
  return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()