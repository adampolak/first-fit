# EVOLVE-BLOCK-START

def construct_intervals():
    """
    Build a sequence of open intervals that forces FirstFit
    to use an exponentially large number of colors while
    the maximum clique size grows only linearly.

    Returns:
      intervals: list of (l, r) tuples in arrival order.
    """
    # recursion depth: increase to raise FF/opt ratio
    k = 5

    # initial seed: one small interval
    T = [(0.0, 1.0)]

    # fixed copy offsets and blocker ranges per iteration
    copy_starts = (2, 6, 10, 14)
    blocker_ranges = ((1, 5), (12, 16), (4, 9), (8, 13))

    for _ in range(k):
        lo = min(l for l, r in T)
        hi = max(r for l, r in T)
        delta = hi - lo

        newT = []
        # replicate T at four scaled positions
        for cs in copy_starts:
            off = delta * cs - lo
            for l, r in T:
                newT.append((l + off, r + off))
        # insert four global blockers that overlap all four copies
        for a, b in blocker_ranges:
            newT.append((delta * a, delta * b))

        T = newT

    return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()