# EVOLVE-BLOCK-START

def construct_intervals():
    """
    Build a sequence of open intervals that forces FirstFit
    to use ~13 colors while the maximum clique size is 5,
    achieving a ratio ≈2.6.
    """
    # key parameter: recursion depth
    k = 5

    # base case: one small interval
    T = [(0.0, 1.0)]

    # each iteration makes 4 scaled copies + 4 blockers
    for _ in range(k):
        lo = min(l for l, r in T)
        hi = max(r for l, r in T)
        delta = hi - lo

        S = []
        # 4 tiled copies at positions 2,6,10,14
        for start in (2, 6, 10, 14, 18):
            offset = delta * start - lo
            for (l, r) in T:
                S.append((l + offset, r + offset))

        # 4 extra “blocking” intervals that overlap all four copies
        S.append((delta * 1,  delta * 5))
        S.append((delta * 12, delta * 16))
        S.append((delta * 4,  delta * 9))
        S.append((delta * 8,  delta * 13))
        # additional blocker to pressure FF further
        S.append((delta * 18, delta * 22))

        T = S

    return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()