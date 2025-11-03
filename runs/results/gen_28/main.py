# EVOLVE-BLOCK-START

def construct_intervals(iterations: int = 4):
    """
    Build a sequence of open intervals that forces FirstFit
    to use a large number of colors while keeping omega small.
    This version parameterizes recursion depth and diversifies patterns
    across iterations to avoid regular structure.
    """
    # base case: a minimal seed, still two? We'll keep single
    T = [(0.0, 1.0)]

    # pattern libraries for diversification
    start_patterns = [
        (2, 6, 10, 14),
        (1, 5, 9, 13),
        (3, 7, 11, 15),
        (0, 4, 8, 12)
    ]

    blocker_templates = [
        [(1, 5), (12, 16), (4, 9), (8, 13)],
        [(0, 4), (11, 15), (3, 8), (7, 12)],
        [(1.5, 5.5), (12.5, 16.5), (4.5, 9.5), (8.5, 13.5)],
        [(0.5, 4.5), (11.5, 15.5), (3.5, 8.5), (7.5, 12.5)]
    ]

    for i in range(iterations):
        lo = min(l for l, r in T)
        hi = max(r for l, r in T)
        delta = hi - lo

        S = []
        starts = start_patterns[i % len(start_patterns)]
        blockers = blocker_templates[i % len(blocker_templates)]

        # Place four scaled/translated copies at staggered offsets (with current pattern)
        for start in starts:
            offset = delta * start - lo
            S.extend([(offset + l, offset + r) for l, r in T])

        # Add four connector intervals using the selected template
        for a, b in blockers:
            S.append((delta * a, delta * b))

        T = S

    return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()