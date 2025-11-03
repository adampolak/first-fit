# EVOLVE-BLOCK-START

def construct_intervals(depth=None):
    """
    Adaptive cap construction with cyclic starts, bridges, and scaling.
    Returns a list of open intervals (l,r) in presentation order.
    """
    # 1) depth control
    if depth is None:
        depth = 5

    # 2) multiple initial seeds to diversify early color usage
    T = [
        (0.0, 1.0),
        (2.0, 3.0)
    ]

    # 3) four rotating start sets (offset multipliers)
    START_SETS = [
        [2, 6, 10, 14],
        [3, 7, 11, 15],
        [4, 8, 12, 16],
        [5, 9, 13, 17],
    ]

    # 4) four rotating bridge patterns
    BRIDGE_SETS = [
        [(1,5),(12,16),(4,9),(8,13)],
        [(2,6),(13,17),(5,10),(9,14)],
        [(3,7),(14,18),(6,11),(10,15)],
        [(4,8),(15,19),(7,12),(11,16)],
    ]

    # 5) non-linear gamma schedule to vary delta each level
    SCALING_LIST = [1.0, 1.25, 0.9, 1.15, 1.05]

    for k in range(depth):
        lo = min(l for l, r in T)
        hi = max(r for l, r in T)
        # cycle gamma perturbations
        gamma = SCALING_LIST[k % len(SCALING_LIST)]
        delta = (hi - lo) * gamma

        S = []
        # cycle start offsets
        start_list = START_SETS[k % len(START_SETS)]
        for s in start_list:
            for (l, r) in T:
                S.append((delta * s + (l - lo), delta * s + (r - lo)))
        # cycle bridges to interlink copies
        bridge_list = BRIDGE_SETS[k % len(BRIDGE_SETS)]
        for (a, b) in bridge_list:
            S.append((delta * a, delta * b))

        T = S

    return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()