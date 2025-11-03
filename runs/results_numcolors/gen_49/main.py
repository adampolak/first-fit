# EVOLVE-BLOCK-START

def construct_intervals(rounds=3, seed_lo=0.0, seed_scale=1.0):
  """
  Deterministic rotating multi-round wave construction with sparse caps and a template bank.
  - rounds: number of recursive expansion rounds
  - seed_lo / seed_scale: deterministic knobs for initial placement and scaling
  Returns:
    intervals: list of tuples (l, r) representing open intervals
  """
  # Start with a compact multi-interval core (promotes early coupling)
  T = [
      (seed_lo, seed_lo + 1.0 * seed_scale),
      (seed_lo + 2.0 * seed_scale, seed_lo + 3.0 * seed_scale),
      (seed_lo + 4.0 * seed_scale, seed_lo + 5.0 * seed_scale),
      (seed_lo + 6.0 * seed_scale, seed_lo + 7.0 * seed_scale),
  ]

  # Four-interval gadget templates (rotate across rounds)
  templates = [
    [(1, 5), (12, 16), (4, 9), (8, 13)],              # A
    [(0.5, 4.5), (11, 15), (3.5, 8.5), (7, 12)],      # B
    [(1, 4), (6, 9), (3, 7), (9, 13)],                # C
    [(2, 6), (7, 11), (0, 3), (10, 14)],              # D
  ]

  # Deterministic cycle of start-patterns: compact 4-wave, dense 8-wave, staggered 7-wave
  start_patterns = [
    (2, 6, 10, 14),
    (2, 4, 6, 8, 10, 12, 14, 16),
    (2, 5, 8, 11, 14, 17, 20),
  ]

  for round_idx in range(rounds):
    # Current span
    lo = min(l for l, r in T)
    hi = max(r for l, r in T)
    delta = hi - lo if hi > lo else 1.0

    S = []

    # Pick starts and gadget for this round
    starts = list(start_patterns[round_idx % len(start_patterns)])
    template = templates[round_idx % len(templates)]

    # Translate copies with alternating inner order; add light bridges
    for i, start in enumerate(starts):
      base = T if ((round_idx + i) % 2 == 0) else list(reversed(T))
      S += [(delta * start + l - lo, delta * start + r - lo) for (l, r) in base]

      if i + 1 < len(starts):
        nxt = starts[i + 1]
        # Light bridge couples colors across adjacent blocks (keeps omega small)
        S.append((delta * (start + 0.5), delta * (nxt + 0.5)))

    # Global gadget intervals scaled by delta
    S += [(delta * a, delta * b) for (a, b) in template]

    # Micro-gadget inside the middle block (scaled down to limit omega)
    middle = starts[len(starts) // 2]
    micro_template = templates[(round_idx + 1) % len(templates)]
    S += [(delta * (middle + a / 16.0), delta * (middle + b / 16.0)) for (a, b) in micro_template]

    # Sparse caps placed relative to block indices to push FF without raising omega much
    max_start = max(starts)
    caps_step = max(3, int(max_start // 4))
    for j in range(2, int(max_start), caps_step):
      S.append((delta * (j - 0.5), delta * (j + 1.5)))

    T = S

  return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()