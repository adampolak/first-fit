# EVOLVE-BLOCK-START

from math import gcd

def normalize_intervals(intervals):
    if not intervals:
        return []
    # scale by 2 to clear halves, translate to ≥0, and divide by gcd
    scaled = [(int(round(l*2)), int(round(r*2))) for l, r in intervals]
    minv = min(l for l, _ in scaled)
    shifted = [ (l-minv, r-minv) for l, r in scaled ]
    vals = [abs(v) for pair in shifted for v in pair]
    g = 0
    for v in vals:
        g = gcd(g, v)
    if g > 1:
        shifted = [(l//g, r//g) for l, r in shifted]
    return shifted

def construct_intervals(iterations=5):
    """
    Recursive adversarial construction:
    - First iteration uses 5 copies to diversify
    - Subsequent iterations use 4 copies
    - All copies and blockers are center-anchored
    """
    # base seed
    T = [(0.0, 1.0)]
    # fixed blocker template
    blockers = [(1.0,5.0), (4.0,9.0), (8.0,13.0), (12.0,16.0)]
    for i in range(iterations):
        lo = min(l for l, _ in T)
        hi = max(r for _, r in T)
        delta = hi - lo
        center = (lo + hi)/2.0
        # dynamic offsets: extra branch first iteration
        if i == 0:
            starts = (2,6,10,14,18)
        else:
            starts = (2,6,10,14)
        S = []
        # make copies, center-anchored
        for s in starts:
            off = delta*s - center
            for l,r in T:
                S.append((l+off, r+off))
        # add blockers, center-anchored
        for a,b in blockers:
            S.append((delta*a - center, delta*b - center))
        T = S
    # normalize and return
    return normalize_intervals(T)

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()