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

    # Mild jitter to break symmetry without inflating omega
    jitter_amp = 0.0 if round_idx == 0 else 0.15

    # Build translated copies with alternating inner order to hamper color reuse
    for idx, start in enumerate(starts):
      # deterministic jitter to avoid self-similarity; in [-jitter_amp, jitter_amp)
      if jitter_amp > 0.0:
        jitter = (((round_idx + 1) * 29 + (idx + 3) * 41) % 101) / 101.0
        jitter = jitter * (2 * jitter_amp) - jitter_amp
      else:
        jitter = 0.0
      adj_start = start + jitter

      base = T if (idx % 2 == 0) else list(reversed(T))
      S += [(delta * adj_start + l - lo, delta * adj_start + r - lo) for l, r in base]

      # Consecutive bridges to couple colors locally
      if idx + 1 < len(starts):
        nxt = starts[idx + 1]
        S.append((delta * (adj_start + 0.45), delta * (nxt + 0.55)))

      # Sparse second-next bridges to couple layers without dense stacking
      if (idx + 2 < len(starts)) and (round_idx % 2 == 0) and ((idx // 2) % 2 == 0):
        nxt2 = starts[idx + 2]
        S.append((delta * (adj_start + 0.70), delta * (nxt2 + 0.30)))

    # Add the chosen template intervals scaled by delta (classic 4-gadget)
    S += [(delta * a, delta * b) for (a, b) in template]

    # Single ladder cap across 3-block span to couple distant colors
    if len(starts) >= 4:
      a = starts[0]
      d = starts[3]
      S.append((delta * (a - 0.20), delta * (d + 0.20)))

    # Micro-caps near a subset of starts (every other) to push FF locally
    cap_len = max(0.6, 0.08 * delta)
    cap_rad = max(0.15, 0.02 * delta)
    for i, start in enumerate(starts):
      if i % 2 == 0:
        mid = delta * start
        S.append((mid - cap_rad, mid + cap_len))

    # Add a pair of long spines only in the middle round to keep omega controlled
    if rounds >= 2 and round_idx == 1:
      left = min(starts) - 0.25
      right = max(starts) + 0.25
      S.append((delta * left, delta * right))
      S.append((delta * (left + 0.20), delta * (right + 0.40)))

    # A few sparse caps across the block range (very low density)
    max_start = max(starts)
    caps_step = max(4, int(max(4, max_start) // 3))
    for j in range(2, max_start, caps_step):
      S.append((delta * (j - 0.40), delta * (j + 1.40)))

    T = S

  return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()