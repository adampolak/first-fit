# EVOLVE-BLOCK-START

def construct_intervals():
    """
    Hybrid two-phase recursive interval construction to force FirstFit high color use.
    """
    # Phase 1 seed: a spine of 4 non-overlapping unit intervals
    T = [(0.0, 1.0), (2.0, 3.0), (4.0, 5.0), (6.0, 7.0)]

    # Patterns for 4-branch starts and a small bank of 4-interval connector gadgets
    start_patterns = [
        (2, 6, 10, 14),
        (1, 5, 9, 13),
        (3, 7, 11, 15),
        (2, 4, 8, 12),
    ]
    template_bank = [
        ((1,5),  (12,16), (4,9),  (8,13)),
        ((0,4),  (5,9),   (2,6),  (7,11)),
        ((1,4),  (6,9),   (3,7),  (9,13)),
        ((2,6),  (7,11),  (0,3),  (10,14)),
    ]

    # Phase 1: major 4-branch expansions
    depth1 = 5
    for i in range(depth1):
        lo = min(l for l, r in T)
        hi = max(r for l, r in T)
        span = hi - lo

        # choose and possibly reverse the start pattern
        pat = start_patterns[i % len(start_patterns)]
        starts = pat if (i % 2 == 0) else pat[::-1]

        S = []
        # create 4 scaled copies, alternating clone order
        for idx, st in enumerate(starts):
            base = span * st - lo
            seq = T if (idx % 2 == 0) else list(reversed(T))
            for (l, r) in seq:
                S.append((l + base, r + base))

        # add exactly 4 connectors from the template bank
        for (a, b) in template_bank[i % len(template_bank)]:
            S.append((span * a, span * b))

        T = S

    # Phase 2: mini-wave expansions at one-eighth scale
    depth2 = 3
    mini_starts = [1, 3, 5, 7]  # fractional offsets over an 8-unit grid
    mini_connectors = [
        ((2/8.0, 6/8.0), (4/8.0, 9/8.0)),
        ((1/8.0, 5/8.0), (3/8.0, 7/8.0)),
    ]

    for j in range(depth2):
        lo2 = min(l for l, r in T)
        hi2 = max(r for l, r in T)
        span2 = hi2 - lo2

        S2 = []
        # inject short‐interval waves
        for st in mini_starts:
            base2 = span2 * (st / 8.0) - lo2
            for (l, r) in T:
                S2.append((l + base2, r + base2))

        # two short cross-connectors per round
        (a1, b1), (a2, b2) = mini_connectors[j % len(mini_connectors)]
        S2.append((span2 * a1 - lo2, span2 * b1 - lo2))
        S2.append((span2 * a2 - lo2, span2 * b2 - lo2))

        T = S2

    # Return the final presentation order
    return [(float(l), float(r)) for (l, r) in T]

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()