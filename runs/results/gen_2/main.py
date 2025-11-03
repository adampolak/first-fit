# EVOLVE-BLOCK-START

def construct_intervals():
    """
    Build a sequence of open intervals that forces FirstFit
    to use many colors while keeping the maximum clique size small.

    Improvements:
    - Increased recursion depth (k=5) to amplify the number of colors FF uses.
    - Normalize final coordinates to keep endpoint magnitudes moderate
      while preserving all intersection relationships.
    """
    # key parameter: recursion depth (increased to intensify the construction)
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
        for start in (2, 6, 10, 14):
            offset = delta * start - lo
            for (l, r) in T:
                S.append((l + offset, r + offset))

        # 4 extra blocking intervals that overlap all four copies
        S.append((delta * 1,  delta * 5))
        S.append((delta * 12, delta * 16))
        S.append((delta * 4,  delta * 9))
        S.append((delta * 8,  delta * 13))

        T = S

    # Normalize coordinates to keep numbers moderate while preserving intersections.
    lo = min(l for l, r in T)
    hi = max(r for l, r in T)
    span = hi - lo if (hi - lo) > 0 else 1.0
    scale = 100.0  # map final span to [0, scale] to keep endpoints tidy
    T = [((l - lo) / span * scale, (r - lo) / span * scale) for (l, r) in T]

    return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()