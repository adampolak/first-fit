# EVOLVE-BLOCK-START

def construct_intervals(rounds=3):
  """
  Construct a sequence of open intervals presented to FirstFit that aims to
  maximize FF colors used relative to the optimal (omega), using a deterministic
  rotating multi-scale expansion.

  Inputs:
    - rounds: number of main-round expansions (default 3)

  Returns:
    - list of (l, r) tuples representing open intervals
  """

  # --------------- Tunable parameters (kept deterministic) ----------------
  # Bridge: short connectors to couple colors across adjacent blocks
  bridge_len_frac = 0.18     # length of bridge relative to delta
  # Caps per round: tiny intervals that discourage color reuse
  caps_per_round = 2
  # Fine phase: small-scale expansion after main rounds
  fine_phase_rounds = 1
  fine_scale_div = 2.5       # delta2 = delta / fine_scale_div
  # Micro-gadget scale inside a selected block
  micro_scale = 1.0 / 16.0
  # -----------------------------------------------------------------------

  # Seed with multiple disjoint unit intervals to promote early diversity
  T = [(0, 1), (2, 3), (4, 5), (6, 7)]

  # Four-interval gadget templates (A, B, C, D)
  templates = [
    [(1, 5), (12, 16), (4, 9), (8, 13)],        # A
    [(0.5, 4.5), (11, 15), (3.5, 8.5), (7, 12)],# B
    [(1, 4), (6, 9), (3, 7), (9, 13)],          # C
    [(2, 6), (7, 11), (0, 3), (10, 14)],        # D
  ]

  # Rotating start patterns (Recommendation 1)
  rotating_starts = [
    (2, 6, 10, 14),
    (1, 5, 9, 13),
    (3, 7, 11, 15),
    (2, 4, 8, 12),
  ]

  # Additional patterns to occasionally densify structure without blowup
  extra_patterns = [
    (2, 4, 6, 8, 10, 12),   # moderate 6-wave
    (2, 5, 8, 11, 14),      # staggered 5-wave
  ]

  def span_of(intervals):
    lo = min(l for l, r in intervals)
    hi = max(r for l, r in intervals)
    return lo, hi, max(hi - lo, 1.0)

  # Main rounds
  for rd in range(max(1, int(rounds))):
    lo, hi, delta = span_of(T)
    S = []

    # Choose starts: rotate primary, then inject an extra pattern every 2 rounds
    base_starts = list(rotating_starts[rd % len(rotating_starts)])
    if rd % 2 == 1:
      base_starts = list(extra_patterns[(rd // 2) % len(extra_patterns)])

    # Select a template (cycle)
    template = templates[rd % len(templates)]

    # Decide block order to break inner symmetry
    if rd % 3 == 0:
      block_order = base_starts
    elif rd % 3 == 1:
      block_order = list(reversed(base_starts))
    else:
      block_order = base_starts[1:] + base_starts[:1]

    # Clone translated copies with alternating inner order
    for bi, st in enumerate(block_order):
      base = T if ((rd + bi) % 2 == 0) else list(reversed(T))
      base_off = delta * st - lo
      for (l, r) in base:
        S.append((l + base_off, r + base_off))

      # Short bridge to next block
      if bi + 1 < len(block_order):
        nxt = block_order[bi + 1]
        b_lo = delta * (st + 0.35)
        b_hi = b_lo + max(delta * bridge_len_frac, 0.5)
        S.append((b_lo, b_hi))

    # Place gadget template scaled by delta
    S += [(delta * a, delta * b) for (a, b) in template]

    # Micro-gadget injected inside a middle block
    mid = block_order[len(block_order) // 2]
    micro_t = templates[(rd + 1) % len(templates)]
    for (a, b) in micro_t:
      S.append((delta * (mid + micro_scale * a), delta * (mid + micro_scale * b)))

    # Sparse towers: couple across two adjacent blocks without raising omega too much
    for bi in range(len(block_order) - 1):
      st = block_order[bi]
      # three staggered layers spanning < 2 blocks
      for layer in range(3):
        off = 0.18 * layer
        t_lo = delta * (st + 0.22 + off)
        t_hi = delta * (st + 1.92 + off)
        if t_hi > t_lo:
          S.append((t_lo, t_hi))

    # Tiny caps to nudge FF
    max_start = max(block_order) if block_order else 2
    positions = [int(max_start * (i + 1) / (caps_per_round + 1)) for i in range(caps_per_round)]
    for pos in positions:
      c_lo = delta * (pos + 0.12)
      c_hi = c_lo + max(0.07 * delta, 0.4)
      S.append((c_lo, c_hi))

    T = S

  # Fine-scale second phase (Recommendation 3)
  if fine_phase_rounds > 0 and len(T) > 0:
    lo, hi, delta = span_of(T)
    delta2 = max(1.0, delta / fine_scale_div)

    # Use a compact 3-wave fine pattern to avoid blowup
    fine_patterns = [(1, 4, 7), (2, 5, 8)]
    for fp in range(fine_phase_rounds):
      starts = list(fine_patterns[fp % len(fine_patterns)])
      S = []
      for bi, st in enumerate(starts):
        # scale down T into the fine blocks (keep topology, reduce span)
        base_off = delta2 * st - lo
        # linear map: x -> (x - lo) * (delta2 / delta) + base_off
        scale = (delta2 / delta)
        base = T if ((fp + bi) % 2 == 0) else list(reversed(T))
        for (l, r) in base:
          S.append(( (l - lo) * scale + base_off, (r - lo) * scale + base_off ))
        # tiny bridge inside fine scale
        if bi + 1 < len(starts):
          b_lo = delta2 * (st + 0.33)
          b_hi = b_lo + max(0.05 * delta2, 0.25)
          S.append((b_lo, b_hi))

      # cap-coupler tying ends within fine scale
      S.append((delta2 * (starts[0] + 0.4), delta2 * (starts[-1] - 0.4)))
      T = S

  return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()