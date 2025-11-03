# EVOLVE-BLOCK-START

def construct_intervals(rounds=3, seed_lo=0.0, seed_scale=1.0):
  """
  Deterministic two-phase construction that pressures FirstFit with layered,
  staggered blocks while heuristically keeping omega small.

  Returns:
    intervals: list of (l, r) tuples (open intervals) in presentation order.
  """

  # -------------------- Configuration and Banks --------------------
  # Hard cap on total intervals; leave headroom for the micro-phase.
  MAX_INTERVALS = 9200
  MICRO_HEADROOM = 350

  # 8-round start-pattern schedule (diversifies positions deterministically)
  start_bank = [
      (2, 6, 10, 14),  # A
      (1, 5, 9, 13),   # B
      (3, 7, 11, 15),  # C
      (2, 4, 8, 12),   # D
      (4, 8, 12, 16),  # E
      (5, 9, 13, 17),  # F
      (3, 5, 9, 13),   # G
      (4, 6, 10, 14),  # H (repeat-like, helps cadence)
  ]

  # 4-interval delta-block connector templates (Figure-4 style variants)
  template_bank = [
      ((1.0, 5.0),  (12.0, 16.0), (4.0, 9.0),  (8.0, 13.0)),  # T0
      ((0.5, 4.5),  (11.0, 15.0), (3.5, 8.5),  (7.0, 12.0)),  # T1
      ((1.0, 4.0),  (6.0, 9.0),   (3.0, 7.0),  (9.0, 13.0)),  # T2
      ((2.0, 6.0),  (7.0, 11.0),  (0.0, 3.0),  (10.0, 14.0)), # T3
  ]

  # Bridge placement (relative to a block start) and bridge min-absolute length
  BRIDGE_FRAC_INNER = 0.42
  BRIDGE_MIN_LEN = 0.18

  # Micro-cap placements for final pressure
  MICRO_CAP_FRACTIONS = (0.18, 0.35, 0.50, 0.65, 0.82)
  MICRO_CAP_WIDTH_FRAC = 0.05  # keep micro caps short to avoid raising omega too much

  # -------------------- Helper Functions --------------------
  def bounds(intervals):
    lo = min(l for l, _ in intervals)
    hi = max(r for _, r in intervals)
    span = hi - lo if hi > lo else 1.0
    return lo, hi, span

  def translate_block(intervals, shift):
    return [(l + shift, r + shift) for (l, r) in intervals]

  def interleave_blocks(blocks, round_idx):
    if not blocks:
      return []
    maxlen = max(len(b) for b in blocks)
    order = list(range(len(blocks)))
    # Change order every two rounds deterministically
    if (round_idx // 2) % 2 == 1:
      order = order[::-1]
    out = []
    for i in range(maxlen):
      for idx in order:
        b = blocks[idx]
        if i < len(b):
          out.append(b[i])
    return out

  def add_template_connectors(span, template, acc):
    # Scale the 4-interval template by current span
    for (a, b) in template:
      acc.append((span * a, span * b))

  def add_kt_connectors(starts, span, acc):
    # Kierstead–Trotter style 4 connectors on start tuple s0..s3 (assumes len >= 4)
    s0, s1, s2, s3 = starts[:4]
    acc.append((span * (s0 - 1), span * (s1 - 1)))  # left cap
    acc.append((span * (s2 + 2), span * (s3 + 2)))  # right cap
    acc.append((span * (s0 + 2), span * (s2 - 1)))  # cross 1
    acc.append((span * (s1 + 2), span * (s3 - 1)))  # cross 2

  def add_short_bridges_between_adjacent(starts, lo, span, acc):
    # Short interior bridges to couple colors across consecutive blocks
    for i in range(len(starts) - 1):
      s = starts[i]
      s_next = starts[i + 1]
      b_lo = span * (s + BRIDGE_FRAC_INNER) - lo
      b_hi = span * (s_next + BRIDGE_FRAC_INNER) - lo
      left, right = (min(b_lo, b_hi), max(b_lo, b_hi))
      # shorten to avoid giant overlaps
      if right - left < BRIDGE_MIN_LEN:
        right = left + BRIDGE_MIN_LEN
      acc.append((left + 0.01 * span, right - 0.01 * span))

  def add_sparse_caps(starts, span, acc):
    # A very light cap family: a few caps spaced across the start lattice
    # Keep very sparse to avoid increasing clique too much.
    if not starts:
      return
    max_start = max(starts)
    step = max(3, int(max_start // 4))
    for j in range(2, max_start, step):
      a = span * (j - 0.45)
      b = span * (j + 0.55)
      acc.append((a + 0.02 * span, b - 0.02 * span))

  def add_micro_caps_global(intervals, acc):
    lo, hi, span = bounds(intervals)
    w = max(MICRO_CAP_WIDTH_FRAC * span, 0.3)
    for frac in MICRO_CAP_FRACTIONS:
      c = lo + frac * span
      acc.append((c - 0.5 * w, c + 0.5 * w))

  def do_round(T, round_idx, starts, template):
    lo, hi, span = bounds(T)
    blocks = []
    # Build translated blocks; alternate reversal parity to hamper reuse
    for bi, s in enumerate(starts):
      base = T[::-1] if ((round_idx + bi) % 2 == 1) else T
      shift = span * s - lo
      blocks.append(translate_block(base, shift))
    # Interleave in round-robin to mix active colors
    S = interleave_blocks(blocks, round_idx)
    # Deterministic connectors (K-T) and a scaled template
    add_kt_connectors(starts, span, S)
    add_template_connectors(span, template, S)
    # Short bridges between adjacent blocks (interior)
    add_short_bridges_between_adjacent(starts, lo, span, S)
    return S

  def normalize_nonnegative(intervals):
    if not intervals:
      return intervals
    min_l = min(l for l, _ in intervals)
    if min_l < 0:
      return [(l - min_l, r - min_l) for (l, r) in intervals]
    return intervals

  # -------------------- Phase 0: Seed --------------------
  # Seed with a 4-interval "spine" to raise early color variety without large omega.
  T = [
      (seed_lo + 0.0 * seed_scale, seed_lo + 1.0 * seed_scale),
      (seed_lo + 2.0 * seed_scale, seed_lo + 3.0 * seed_scale),
      (seed_lo + 4.0 * seed_scale, seed_lo + 5.0 * seed_scale),
      (seed_lo + 6.0 * seed_scale, seed_lo + 7.0 * seed_scale),
  ]

  # -------------------- Phase 1: Coarse multiscale expansion --------------------
  # Adaptive number of rounds: aim near MAX_INTERVALS - MICRO_HEADROOM
  # With 4x replication per round and a few connectors, 5 rounds from 4 seeds is safe.
  target_limit = MAX_INTERVALS - MICRO_HEADROOM
  coarse_rounds = max(5, int(rounds) + 2)  # deterministic; default rounds=3 -> 5 coarse rounds

  for ridx in range(coarse_rounds):
    starts = start_bank[ridx % len(start_bank)]
    template = template_bank[ridx % len(template_bank)]
    T = do_round(T, ridx, starts, template)
    # Light sparse caps after each round (kept sparse)
    lo, hi, span = bounds(T)
    add_sparse_caps(starts, span, T)
    # Safety guard
    if len(T) >= target_limit:
      # trim deterministically from the front (keep recent structure for late pressure)
      T = T[-target_limit:]
      break

  # -------------------- Phase 2: Light micro delta coupling --------------------
  # Add a thin second phase that adds a few scaled template gadgets and bridges,
  # not full clones of T (to avoid blowing up size).
  lo, hi, span = bounds(T)
  # Place 3 anchor starts inside the global span for micro gadgets
  micro_starts = (1.5, 4.5, 7.5)
  micro_scale = 0.18  # small scale to keep omega modest

  S2 = list(T)  # start from current T; append micro structures
  for k, s in enumerate(micro_starts):
    # Pick a template and place it near start s at micro scale
    template = template_bank[(k + coarse_rounds) % len(template_bank)]
    base = span * s
    for (a, b) in template:
      S2.append((base + span * micro_scale * a, base + span * micro_scale * b))
    # short bridge within this micro site
    S2.append((base + 0.25 * span * micro_scale, base + 0.75 * span * micro_scale))

  # Lastly, sprinkle micro caps across the entire geometry (late arrivals)
  add_micro_caps_global(T, S2)
  T = S2

  # Final normalization and safety trim
  T = normalize_nonnegative(T)
  if len(T) > MAX_INTERVALS:
    T = T[-MAX_INTERVALS:]

  return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()