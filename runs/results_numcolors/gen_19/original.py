# EVOLVE-BLOCK-START

def construct_intervals():
    """
    Construct a deterministic sequence of open intervals in presentation order.
    This version keeps the same multiplicative 4-branch recursion (so the
    total count stays comparable to the strong prior examples) but
    diversifies arrival order, clone order, and connector templates across
    iterations to make FirstFit reuse less effective.

    Returns:
        intervals: list of tuples (l, r) representing open intervals.
    """
    # Target recursion depth (6 is a practical sweet spot: ~9556 intervals for
    # the classic 4-branch + 4-connectors growth). Keep fixed for reproducibility.
    depth = 6

    # Base seed
    T = [(0.0, 1.0)]

    # A cycle of 4 start-patterns (each has 4 block starts -> branching=4)
    start_patterns = [
        (2, 6, 10, 14),
        (1, 5, 9, 13),
        (3, 7, 11, 15),
        (2, 4, 8, 12),
    ]

    # A small bank of 4-connector templates (each template adds exactly 4 intervals)
    # We pick one per iteration to vary coupling without increasing per-iteration cost.
    template_bank = [
        ((1, 5), (12, 16), (4, 9), (8, 13)),   # classical gadget
        ((0, 4), (5, 9), (2, 6), (7, 11)),
        ((1, 4), (6, 9), (3, 7), (9, 13)),
        ((2, 6), (7, 11), (0, 3), (10, 14)),
    ]

    for it in range(depth):
        lo = min(l for l, r in T)
        hi = max(r for l, r in T)
        span = hi - lo

        # pick a start pattern from the cycle
        starts = list(start_patterns[it % len(start_patterns)])

        # vary the order in which blocks are appended to S; this changes arrival order
        if it % 3 == 0:
            block_order = starts
        elif it % 3 == 1:
            block_order = list(reversed(starts))
        else:
            # rotate by 1 to alter interleaving
            block_order = starts[1:] + starts[:1]

        S = []
        for idx, st in enumerate(block_order):
            # alternate clone insertion order to break symmetry:
            # sometimes append copies of T in forward order, sometimes reversed.
            if (it + idx) % 2 == 0:
                clone_order = T
            else:
                clone_order = list(reversed(T))

            base_off = span * st - lo
            for (l, r) in clone_order:
                S.append((l + base_off, r + base_off))

        # append exactly 4 connectors chosen from the template bank (keeps per-iteration
        # additive cost constant so total size remains in the same growth family)
        template = template_bank[it % len(template_bank)]
        for (a, b) in template:
            S.append((span * a, span * b))

        # update for next iteration
        T = S

    return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()