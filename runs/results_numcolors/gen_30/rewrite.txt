# EVOLVE-BLOCK-START

def construct_intervals(rounds=3):
  """
  Construct a sequence of intervals of real line,
  in the order in which they are presented to FirstFit,
  so that it tends to maximize the number of colors used by FirstFit
  divided by the maximum number of intervals that cover a single point
  through a recursive wave construction.

  Improvements blended via crossover:
  - multi-unit seed (keeps omega small, encourages coupling)
  - cycled start patterns (varied offsets across rounds)
  - alternating inner order (hampers color reuse)
  - deterministic bridging (consecutive and second-next) to couple colors
  - per-round four-interval gadget templates
  - micro-caps and sparse long caps to build towers without large omega
  - moderate growth control to keep size manageable
  """
  # Seed with several disjoint unit intervals (open intervals)
  T = [(0, 1), (2, 3), (4, 5), (6, 7)]

  # Template bank (scaled each round)
  templates = [
    [(1, 5), (12, 16), (4, 9), (8, 13)],          # A: classic
    [(0.5, 4.5), (11, 15), (3.5, 8.5), (7, 12)],  # B: shifted
    [(1, 4), (6, 9), (3, 7), (9, 13)],            # C: tighter overlaps
    [(2, 6), (7, 11), (0, 3), (10, 14)],          # D: staggered
  ]

  # Start patterns: cycle deterministically
  start_patterns = [
    (2, 6, 10, 14),                    # compact 4-wave
    (2, 4, 6, 8, 10, 12, 14, 16),      # denser 8-wave
    (2, 5, 8, 11, 14, 17, 20),         # staggered 7-wave
  ]

  for round_idx in range(rounds):
    lo = min(l for l, r in T)
    hi = max(r for l, r in T)
    delta = hi - lo if hi > lo else 1.0

    starts = list(start_patterns[round_idx % len(start_patterns)])
    template = templates[round_idx % len(templates)]
    S = []

    # Clone T into translated blocks; alternate order; add bridges
    for idx, start in enumerate(starts):
      base = T if (idx % 2 == 0) else list(reversed(T))
      shift = delta * start
      S += [(shift + (l - lo), shift + (r - lo)) for (l, r) in base]

      # Bridge to next block to couple colors locally
      if idx + 1 < len(starts):
        nxt = starts[idx + 1]
        S.append((delta * (start + 0.5), delta * (nxt + 0.5)))

      # Optional bridge to second-next to couple across layers
      if idx + 2 < len(starts) and (idx % 2 == 0):
        nxt2 = starts[idx + 2]
        S.append((delta * (start + 0.75), delta * (nxt2 + 0.25)))

    # Add the chosen template scaled by delta
    S += [(delta * a, delta * b) for (a, b) in template]

    # Micro-caps near each start to press FirstFit locally without raising omega much
    cap_radius = max(0.2, 0.03 * delta)
    cap_length = max(0.8, 0.12 * delta)
    for start in starts:
      mid = delta * start
      # two asymmetric micro-caps per start
      S.append((mid - cap_radius, mid + cap_length))
      S.append((mid - 0.5 * cap_radius, mid + 0.8 * cap_length))

    # Sparse long caps across blocks to build towers/coupling while keeping density low
    max_start = max(starts)
    step = max(3, max_start // 4)
    for j in range(2, max_start, step):
      S.append((delta * (j - 1.0), delta * (j + 2.0)))

    # Light "tower caps": connect every other block (avoid dense stacking)
    for i in range(0, len(starts) - 2, 2):
      a, c = starts[i], starts[i + 2]
      S.append((delta * (a - 0.25), delta * (c + 0.25)))

    T = S

  return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()