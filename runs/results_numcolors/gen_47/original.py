# EVOLVE-BLOCK-START

def construct_intervals(rounds=3, seed_lo=0.0, seed_scale=1.0):
  """
  Deterministic, adaptive multi-round wave construction.
  - rounds: number of recursive expansion rounds
  - seed_lo / seed_scale: deterministic knobs for initial placement and scaling
  Returns:
    intervals: list of intervals (l, r) as open intervals
  """
  # Seed with several disjoint unit intervals to promote early overlap coupling
  T = [(seed_lo, seed_lo + 1.0 * seed_scale),
       (seed_lo + 2.0 * seed_scale, seed_lo + 3.0 * seed_scale),
       (seed_lo + 4.0 * seed_scale, seed_lo + 5.0 * seed_scale),
       (seed_lo + 6.0 * seed_scale, seed_lo + 7.0 * seed_scale)]

  # A small bank of four-interval gadgets (templates)
  templates = [
    # Template A: original pattern
    [(1, 5), (12, 16), (4, 9), (8, 13)],
    # Template B: shifted
    [(0.5, 4.5), (11, 15), (3.5, 8.5), (7, 12)],
    # Template C: tighter internal overlaps
    [(1, 4), (6, 9), (3, 7), (9, 13)],
    # Template D: staggered caps
    [(2, 6), (7, 11), (0, 3), (10, 14)]
  ]

  # Deterministic cycle of start-patterns (each is a tuple of starts)
  start_patterns = [
    (2, 6, 10, 14),                 # compact 4-wave
    (2, 4, 6, 8, 10, 12, 14, 16),   # dense 8-wave
    (2, 5, 8, 11, 14, 17, 20),      # staggered 7-wave
  ]

  for round_idx in range(rounds):
    lo = min(l for l, r in T)
    hi = max(r for l, r in T)
    delta = hi - lo
    if delta <= 0:
      delta = 1.0  # safety

    S = []

    # pick deterministic pattern/template for this round
    starts = list(start_patterns[round_idx % len(start_patterns)])
    template = templates[round_idx % len(templates)]

    # Build translated copies with alternating inner order to hamper color reuse
    for idx, start in enumerate(starts):
      base = T if (idx % 2 == 0) else list(reversed(T))
      S += [(delta * start + l - lo, delta * start + r - lo) for l, r in base]

      # light bridging between consecutive blocks to couple colors
      if idx + 1 < len(starts):
        nxt = starts[idx + 1]
        S.append((delta * (start + 0.5), delta * (nxt + 0.5)))

    # Add the chosen template intervals scaled by delta
    S += [(delta * a, delta * b) for (a, b) in template]

    # Add sparse caps to push color usage without large omega
    max_start = max(starts)
    caps_step = max(3, int(max(3, max_start) // 4))
    for j in range(2, max_start, caps_step):
      # small caps around positions to lightly pressure color usage
      S.append((delta * (j - 0.5), delta * (j + 1.0 + 0.5)))

    T = S

  return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()