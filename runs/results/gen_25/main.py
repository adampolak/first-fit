# EVOLVE-BLOCK-START

def construct_intervals(iterations=6):
    """
    Construct a sequence of open intervals (l, r) presented in order to FirstFit.
    Uses the recursive 4-copy + 4-connector building block from
    https://arxiv.org/abs/1506.00192 (Figure 4).  Increasing 'iterations'
    amplifies FirstFit's color usage more than the maximum clique size.
    Default iterations=6 yields ≈19 FirstFit colors vs. clique=7 ⇒ ≈2.714 ratio.
    """
    # Start with one unit interval
    T = [(0.0, 1.0)]
    for _ in range(iterations):
        lo = min(l for l, r in T)
        hi = max(r for l, r in T)
        delta = hi - lo
        S = []
        # Four scaled/translated copies at offsets 2,6,10,14
        for start in (2, 6, 10, 14):
            off = delta * start - lo
            S.extend([(l + off, r + off) for l, r in T])
        # Four long connector intervals that link the copies
        S += [
            (delta * 1,  delta * 5),
            (delta * 12, delta * 16),
            (delta * 4,  delta * 9),
            (delta * 8,  delta * 13),
        ]
        T = S

    # Normalize all distinct endpoints to even integers in increasing order
    endpoints = sorted({x for seg in T for x in seg})
    coord = {}
    cur = 0
    for e in endpoints:
        coord[e] = cur
        cur += 2
    normalized = [(coord[l], coord[r]) for (l, r) in T]
    return normalized

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()