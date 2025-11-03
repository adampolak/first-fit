# EVOLVE-BLOCK-START

def construct_intervals():
  """
  Construct a sequence of intervals on the real line,
  presented to FirstFit in this order, attempting to maximize
  FF colors divided by the clique number (omega).
  Returns a list of (l, r) open intervals with integer endpoints.
  """

  # Hard capacity to keep total intervals < 10000
  CAP = 9800

  # Global deterministic parameters
  SPINE_TEMPLATES = [
    (2, 6, 10, 14),  # classic KT
    (1, 5, 9, 13),   # left-shifted
    (3, 7, 11, 15),  # right-shifted
    (4, 8, 12, 16),  # stretched
  ]
  DEFAULT_STARTS = (2, 6, 10, 14)

  # Optional density knob for future experimentation (kept at 1 here).
  K_DENSITY = 1  # can be toggled to 2 if evaluated safe; left at 1 for stability

  # Stage 0: seed with one unit interval
  T = [(0, 1)]

  # -------------------------------
  # Shared utilities
  # -------------------------------
  def span_info(seq):
    lo = min(l for l, r in seq)
    hi = max(r for l, r in seq)
    delta = hi - lo
    if delta <= 0:
      delta = 1
    return lo, hi, delta

  def kt_connectors(starts, delta, base_shift=0):
    s0, s1, s2, s3 = starts
    cons = [
      (base_shift + (s0 - 1) * delta, base_shift + (s1 - 1) * delta),  # left cap
      (base_shift + (s2 + 2) * delta, base_shift + (s3 + 2) * delta),  # right cap
      (base_shift + (s0 + 2) * delta, base_shift + (s2 - 1) * delta),  # cross 1
      (base_shift + (s1 + 2) * delta, base_shift + (s3 - 1) * delta),  # cross 2
    ]
    return [(a, b) for (a, b) in cons if b > a]

  def spine_round(current_T, starts, interleave=False, reverse_blocks=False, add_long=False):
    lo, hi, delta = span_info(current_T)
    # Build translated blocks
    blocks = []
    for s in starts:
      base = s * delta - lo
      block = [(l + base, r + base) for (l, r) in current_T]
      blocks.append(block[::-1] if reverse_blocks else block)

    # Merge blocks
    S = []
    if interleave:
      maxlen = max(len(b) for b in blocks)
      order = list(range(len(blocks)))
      for i in range(maxlen):
        for idx in order:
          if i < len(blocks[idx]):
            S.append(blocks[idx][i])
    else:
      for blk in blocks:
        S.extend(blk)

    # Add classic connectors
    S.extend(kt_connectors(starts, delta, base_shift=0))
    # Optionally add a longer connector to couple extremes
    if add_long:
      s0, s3 = starts[0], starts[3]
      long_c = (int((s0 + 3) * delta), int((s3 + 3) * delta))
      if long_c[1] > long_c[0]:
        S.append(long_c)
    return S

  def predicted_next_size(sz):
    # KT growth per round: 4*size + 4
    return 4 * sz + 4

  def thin_seed(current_T, max_seed):
    n = len(current_T)
    if n == 0 or max_seed <= 0:
      return []
    step = max(1, n // max_seed)
    return current_T[::step][:max_seed]

  def enforce_budget(seq, budget):
    if budget <= 0:
      return []
    return seq[:budget] if len(seq) > budget else seq

  # -------------------------------
  # Stage 1: Deterministic KT spine
  # -------------------------------
  # Derive a deterministic per-round interleaving parity with seed = f(round_index)
  # Even rounds: interleave, Odd: sequential; reverse block order on odd rounds for variety.
  # Rotate templates each round to diversify positions while preserving clique control.
  max_rounds = 6  # known to fit comfortably under CAP with some room for micro-phase
  for ridx in range(max_rounds):
    # Predict size; stop if we'd exceed CAP by too much (leave room for micro-phase)
    nxt = predicted_next_size(len(T))
    if nxt > CAP:
      break

    starts = SPINE_TEMPLATES[ridx % len(SPINE_TEMPLATES)]
    # Parity discipline (deterministic, seedless)
    interleave = (ridx % 2 == 0)
    reverse_blocks = (ridx % 2 == 1)
    add_long = (ridx % 2 == 0)
    T = spine_round(T, starts, interleave=interleave, reverse_blocks=reverse_blocks, add_long=add_long)

    # Bail if we unexpectedly hit CAP
    if len(T) >= CAP:
      T = T[:CAP]
      return T

  # If already near CAP, return spine-only baseline
  if len(T) >= CAP - 16:
    return T

  # -------------------------------
  # Stage 1.5: Near-tail sparse caps
  # -------------------------------
  # Insert a few long caps near the tail to raise FF mixing with minimal omega impact.
  def insert_near_tail(seq, extra):
    out = list(seq)
    for i, iv in enumerate(extra):
      pos = len(out) - (i * 2 + 1)
      if pos < 0:
        out.append(iv)
      else:
        out.insert(pos, iv)
    return out

  glo, ghi, G = span_info(T)
  def cap_at(fr_a, fr_b):
    L = glo + max(1, int(round(fr_a * G)))
    R = glo + max(1, int(round(fr_b * G)))
    if R <= L:
      R = L + 1
    return (L, R)

  tail_caps = [cap_at(0.08, 0.60), cap_at(0.25, 0.75), cap_at(0.75, 0.92)]
  room = CAP - len(T)
  if room > 0:
    T = insert_near_tail(T, tail_caps[:room])

  if len(T) >= CAP - 24:
    return T

  # -------------------------------
  # Stage 2: Two guarded micro-phases
  # -------------------------------
  # Primary micro-phase (A): windows with small deterministic shifts
  # Secondary micro-phase (B): distinct alternate windows set, includes a cross4 layer

  def build_micro_round(current_T, budget, iter_id=0, alt=False):
    if not current_T or budget <= 8:
      return []

    glo = min(l for l, r in current_T)
    ghi = max(r for l, r in current_T)
    G = max(1, ghi - glo)

    # Thin, evenly spaced seed — bounded by budget and overall size
    seed_sz = max(8, min(40, len(current_T) // 250 if len(current_T) >= 250 else 16))
    U = thin_seed(current_T, seed_sz)
    if not U:
      return []

    ulo = min(l for l, r in U)

    # Window sets
    if not alt:
      # Primary windows with a tiny, deterministic shift based on iter_id
      shift = (iter_id % 3) * 0.02
      fr_windows = [
        (0.12 + shift, 0.22 + shift),
        (0.35 + shift, 0.45 + shift),
        (0.58 + shift, 0.68 + shift),
        (0.80 + shift, 0.90 + shift),
      ]
    else:
      # Alternate windows (distinct, fixed)
      fr_windows = [(0.05, 0.15), (0.28, 0.38), (0.60, 0.70), (0.82, 0.92)]

    # Clamp windows within (0,1)
    fr_windows = [
      (max(0.02, min(0.90, a)), max(0.10, min(0.98, b))) for (a, b) in fr_windows
    ]

    # Build translated micro-blocks aligned to windows
    blocks = []
    for wi, (fa, fb) in enumerate(fr_windows):
      win_lo = glo + int(round(fa * G))
      base = win_lo - ulo
      block = [(l + base, r + base) for (l, r) in U]
      # Parity-based internal reversal per block for mixing
      tag = iter_id + (1 if alt else 0)
      if ((wi + tag) % 2) == 0:
        block = list(reversed(block))
      blocks.append(block)

    # Interleave micro-blocks with parity discipline
    micro = []
    maxlen = max(len(b) for b in blocks)
    order = list(range(len(blocks)))
    if (iter_id % 2) == 1:
      order.reverse()
    for i in range(maxlen):
      for idx in order:
        blk = blocks[idx]
        if i < len(blk):
          micro.append(blk[i])

    # Deterministic connectors across windows (fractional-span analog of KT caps)
    # Always add the classic four; add an extra long cross4 layer when alt=True.
    micro_connectors = [
      (glo + int(round(0.08 * G)), glo + int(round(0.30 * G))),  # left cap
      (glo + int(round(0.60 * G)), glo + int(round(0.92 * G))),  # right cap
      (glo + int(round(0.26 * G)), glo + int(round(0.56 * G))),  # cross 1
      (glo + int(round(0.44 * G)), glo + int(round(0.78 * G))),  # cross 2
    ]
    if alt:
      # Cross4 layer: longer-range couplers to tie distant windows without dense stacking
      micro_connectors.extend([
        (glo + int(round(0.18 * G)), glo + int(round(0.84 * G))),
        (glo + int(round(0.10 * G)), glo + int(round(0.70 * G))),
        (glo + int(round(0.34 * G)), glo + int(round(0.96 * G))),
      ])

    for a, b in micro_connectors:
      if b > a:
        micro.append((a, b))

    # Enforce budget
    micro = enforce_budget(micro, budget)
    return micro

  # Execute primary micro-phase (A), then secondary (B), both budgeted
  # Phase A: up to two iterations (iter_id = 0, 1)
  for iter_id in range(2):
    room = CAP - len(T)
    if room <= 8:
      break
    microA = build_micro_round(T, budget=room, iter_id=iter_id, alt=False)
    if not microA:
      break
    T.extend(enforce_budget(microA, CAP - len(T)))

  # Phase B: a single alternate-window pass with cross4 layer
  room = CAP - len(T)
  if room > 8:
    microB = build_micro_round(T, budget=room, iter_id=2, alt=True)
    if microB:
      T.extend(enforce_budget(microB, CAP - len(T)))

  # Final cap
  if len(T) > CAP:
    T = T[:CAP]

  return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()