# EVOLVE-BLOCK-START

def construct_intervals(rounds=3, seed_lo=0.0, seed_scale=1.0):
    """
    Deterministic, adaptive multi-round wave construction with delta-scaling,
    template gadgets, patterned starts, and light bridging/caps to increase
    FirstFit colors while keeping the clique number modest.
    Returns:
      intervals: list of (l, r) open intervals in presentation order.
    """
    # Seed: several disjoint unit intervals to keep omega low but enable coupling
    T = [(seed_lo, seed_lo + 1.0 * seed_scale),
         (seed_lo + 2.0 * seed_scale, seed_lo + 3.0 * seed_scale),
         (seed_lo + 4.0 * seed_scale, seed_lo + 5.0 * seed_scale),
         (seed_lo + 6.0 * seed_scale, seed_lo + 7.0 * seed_scale)]

    # Bank of four-interval gadgets (scaled by delta each round)
    templates = [
        [(1, 5), (12, 16), (4, 9), (8, 13)],        # A: classic
        [(0.5, 4.5), (11, 15), (3.5, 8.5), (7, 12)],# B: shifted
        [(1, 4), (6, 9), (3, 7), (9, 13)],          # C: tighter overlaps
        [(2, 6), (7, 11), (0, 3), (10, 14)],        # D: staggered caps
    ]

    # Deterministic cycle of start-patterns
    start_patterns = [
        (2, 6, 10, 14),                 # compact 4-wave
        (2, 4, 6, 8, 10, 12, 14, 16),   # dense 8-wave
        (2, 5, 8, 11, 14, 17, 20),      # staggered 7-wave
    ]

    for round_idx in range(rounds):
        lo = min(l for l, r in T)
        hi = max(r for l, r in T)
        delta = hi - lo if hi > lo else 1.0

        S = []

        starts = list(start_patterns[round_idx % len(start_patterns)])
        template = templates[round_idx % len(templates)]

        # Translated copies with alternating inner order to impede color reuse
        for i, start in enumerate(starts):
            base = T if (i % 2 == 0) else list(reversed(T))
            S += [(delta * start + l - lo, delta * start + r - lo) for (l, r) in base]

            # Light bridges to couple colors across adjacent and occasional second-next blocks
            if i + 1 < len(starts):
                nxt = starts[i + 1]
                S.append((delta * (start + 0.5), delta * (nxt + 0.5)))
            if (round_idx % 2 == 0) and (i + 2 < len(starts)):
                nxt2 = starts[i + 2]
                S.append((delta * (start + 0.75), delta * (nxt2 + 0.25)))

        # Insert the chosen template scaled by delta
        S += [(delta * a, delta * b) for (a, b) in template]

        # Sparse caps to increase FF pressure without large local cliques
        max_start = max(starts)
        step = max(3, int(max_start // 4))
        for j in range(2, max_start, step):
            S.append((delta * (j - 0.5), delta * (j + 1.0 + 0.5)))

        T = S

    return T

def run_experiment(**kwargs):
    # Use 3 rounds by default; can be overridden by kwargs
    return construct_intervals(**({k: v for k, v in kwargs.items() if k in {"rounds", "seed_lo", "seed_scale"}}))

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()