# EVOLVE-BLOCK-START

def construct_intervals(seed_count=1):
  """
  Construct a sequence of intervals on the real line, in the order presented to FirstFit,
  to maximize FF colors divided by the clique number (omega).
  Interface preserved: construct_intervals(seed_count=1) -> list[(l, r)]
  """

  # Hard capacity to respect the 10000 cap and leave room for micro-passes
  CAP = 9800
  BASE_SEED = 137  # deterministic seed anchor (no RNG used)

  # Six spine templates (four-start each). Rotate across rounds to couple scales.
  SPINE_TEMPLATES = [
    (2, 6, 10, 14),  # classic KT
    (1, 5, 9, 13),   # left-shifted
    (3, 7, 11, 15),  # right-shifted
    (4, 8, 12, 16),  # stretched-right
    (2, 5, 9, 14),   # mixed-1
    (1, 6, 10, 15),  # mixed-2
  ]

  # Interleaving schedule for six spine rounds (deterministic pattern)
  INTERLEAVE_FLAGS = [True, False, True, True, False, True]
  REVERSE_FLAGS    = [False, True, False, True, False, True]

  # Windows for the two micro passes (fractional positions over the global span)
  MICRO_WINDOWS_A = [(0.12, 0.22), (0.32, 0.42), (0.52, 0.62), (0.72, 0.82)]
  MICRO_WINDOWS_B = [(0.18, 0.28), (0.40, 0.50), (0.62, 0.72), (0.84, 0.92)]

  # Deterministic long-range connector fractions (inserted right after the spine)
  LONG_CONNECTOR_FRACS = [
    (0.08, 0.30), (0.20, 0.55), (0.38, 0.85),
    (0.12, 0.78), (0.46, 0.92), (0.05, 0.48),
  ]

  # Seed: keep single unit interval; multi-seed increases omega prematurely
  T = [(0, 1)]

  # -------------------------------
  # Utilities (purely deterministic)
  # -------------------------------
  def span_delta(seq):
    lo = min(l for l, r in seq)
    hi = max(r for l, r in seq)
    d = hi - lo
    return lo, hi, (d if d > 0 else 1)

  def translate_block(seq, base):
    # Translate every interval by +base
    return [(l + base, r + base) for (l, r) in seq]

  def interleave_blocks(blocks, reverse_order=False):
    # Interleave four blocks row-wise, optionally reversing block order
    if not blocks:
      return []
    order = list(range(len(blocks)))
    if reverse_order:
      order.reverse()
    maxlen = max(len(b) for b in blocks)
    out = []
    for i in range(maxlen):
      for idx in order:
        blk = blocks[idx]
        if i < len(blk):
          out.append(blk[i])
    return out

  def append_connectors(out, starts, delta):
    # Classic four connectors used in KT-style constructions
    s0, s1, s2, s3 = starts
    out.append(((s0 - 1) * delta, (s1 - 1) * delta))  # left cap
    out.append(((s2 + 2) * delta, (s3 + 2) * delta))  # right cap
    out.append(((s0 + 2) * delta, (s2 - 1) * delta))  # cross 1
    out.append(((s1 + 2) * delta, (s3 - 1) * delta))  # cross 2

  def predict_next_size(sz):
    # Pure KT growth from current size: 4*sz + 4
    return 4 * sz + 4

  def even_sample(seq, k):
    # Take an evenly spaced sample of size up to k (k >= 1)
    if not seq or k <= 0:
      return []
    n = len(seq)
    step = max(1, n // k)
    return [seq[i] for i in range(0, n, step)][:k]

  def frac_interval(lo, delta, a_frac, b_frac):
    L = lo + max(1, int(round(a_frac * delta)))
    R = lo + max(1, int(round(b_frac * delta)))
    if R <= L:
      R = L + 1
    return (L, R)

  # -------------------------
  # Spine builder (six rounds)
  # -------------------------
  def apply_spine_round(current_T, starts, do_interleave, reverse_blocks):
    lo, hi, delta = span_delta(current_T)
    # Build four translated blocks
    blocks = []
    for s in starts:
      base = s * delta - lo
      blocks.append(translate_block(current_T, base))

    # Combine blocks (interleaved or sequential)
    if do_interleave:
      S = interleave_blocks(blocks, reverse_order=reverse_blocks)
    else:
      if reverse_blocks:
        blocks = list(reversed(blocks))
      S = []
      for blk in blocks:
        S.extend(blk)

    # Append KT-style connectors
    append_connectors(S, starts, delta)
    return S

  def build_spine(initial_T, rounds, cap_guard):
    T_local = list(initial_T)
    completed_rounds = 0
    for ridx in range(rounds):
      # Predict growth; keep headroom for tail operations (connectors + micro)
      predicted = predict_next_size(len(T_local))
      if predicted >= cap_guard:
        break
      starts = SPINE_TEMPLATES[ridx % len(SPINE_TEMPLATES)]
      do_inter = INTERLEAVE_FLAGS[ridx % len(INTERLEAVE_FLAGS)]
      rev = REVERSE_FLAGS[ridx % len(REVERSE_FLAGS)]
      T_local = apply_spine_round(T_local, starts, do_inter, rev)
      completed_rounds += 1
      if len(T_local) >= cap_guard:
        T_local = T_local[:cap_guard]
        break
    return T_local, completed_rounds

  # -----------------------------------
  # Deterministic long-range connectors
  # -----------------------------------
  def add_long_range_connectors(seq, budget):
    if budget <= 0:
      return seq
    lo, hi, delta = span_delta(seq)
    conns = []
    for (a, b) in LONG_CONNECTOR_FRACS:
      conns.append(frac_interval(lo, delta, a, b))
    if len(conns) > budget:
      conns = conns[:budget]
    return seq + conns

  # ---------------------------
  # Micro-pass builder (two-pass)
  # ---------------------------
  def build_micro_pass(current_T, budget, windows, pass_id):
    if budget <= 0 or not current_T:
      return []

    glo, ghi, G = span_delta(current_T)
    # Thin seed size is tethered to the spine size; keep small to control omega
    seed_sz = max(8, min(36, len(current_T) // 300))
    U = even_sample(current_T, seed_sz)
    if not U:
      return []

    ulo, _, _ = span_delta(U)

    # Build translated micro-blocks
    blocks = []
    for (fa, fb) in windows:
      win_lo = glo + int(round(fa * G))
      base = win_lo - ulo
      blk = translate_block(U, base)
      # Alternate internal reversal to break symmetry
      if ((int(round(fa * 100)) + pass_id) // 5) % 2 == 1:
        blk = list(reversed(blk))
      blocks.append(blk)

    # Interleave micro-blocks with a deterministic pattern per pass
    interleave_rev = bool(pass_id % 2 == 1)
    micro = interleave_blocks(blocks, reverse_order=interleave_rev)

    # Fractional connectors at micro scale (sparse, omega-friendly)
    micro_connectors = [
      frac_interval(glo, G, 0.08, 0.30),
      frac_interval(glo, G, 0.26, 0.56),
      frac_interval(glo, G, 0.60, 0.92),
    ]
    micro.extend(micro_connectors)

    # Trim to budget
    if len(micro) > budget:
      micro = micro[:budget]
    return micro

  # --------------------------
  # Assembly using new layout
  # --------------------------
  # Phase 1: Build the spine with six rounds and a modest CAP guard (reserve tail headroom)
  # Headroom: reserve ~200 intervals for long-range connectors and micro scaffolding deterministically
  SPINE_HEADROOM = 220
  T, depth_done = build_spine(T, rounds=6, cap_guard=max(1, CAP - SPINE_HEADROOM))

  # Phase 1.5: Deterministic long-range connectors (post-spine)
  room = CAP - len(T)
  if room > 0:
    # Use at most 6 connectors; keep a buffer of 160 for micro phases
    conn_budget = max(0, min(6, room - 160))
    if conn_budget > 0:
      T = add_long_range_connectors(T, conn_budget)

  # Early exit if close to CAP already
  if len(T) >= CAP - 16:
    return T

  # Phase 2: Two-pass micro phases with CAP-aware budgeting
  room = CAP - len(T)
  if room > 0:
    # Split deterministically: 60% to pass A, 40% to pass B (at least a few intervals each)
    budget_A = max(12, (room * 3) // 5)
    budget_B = max(8, room - budget_A)

    micro_A = build_micro_pass(T, budget_A, MICRO_WINDOWS_A, pass_id=0)
    if micro_A:
      take = min(len(micro_A), CAP - len(T))
      T.extend(micro_A[:take])

    # Recompute room before pass B
    room = CAP - len(T)
    if room > 0:
      micro_B = build_micro_pass(T, budget_B, MICRO_WINDOWS_B, pass_id=1)
      if micro_B:
        take = min(len(micro_B), CAP - len(T))
        T.extend(micro_B[:take])

  # Capacity trim (safety)
  if len(T) > CAP:
    T = T[:CAP]

  return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()