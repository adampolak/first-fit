# EVOLVE-BLOCK-START

from math import gcd

def normalize_to_grid(intervals):
    """
    Map unique endpoints to 0,2,4,... and shrink by gcd.
    """
    if not intervals:
        return []
    pts = sorted({x for seg in intervals for x in seg})
    coord = {v: i*2 for i,v in enumerate(pts)}
    norm = [(coord[l], coord[r]) for (l,r) in intervals]
    g = 0
    for l,r in norm:
        g = gcd(g, abs(l)); g = gcd(g, abs(r))
    if g>1:
        norm = [(l//g,r//g) for l,r in norm]
    return norm

def construct_intervals():
    """
    Build a k‐level recursive tiling with cycle patterns and wave augmentation.
    """
    # base seed
    T = [(0.0, 1.0)]
    k = 4
    # cycle among four copy‐offset patterns
    copy_patterns = [
        (2, 6, 10, 14),
        (1, 5, 9, 13),
        (3, 7, 11, 15),
        (0, 4, 8, 12),
    ]
    # cycle among three blocker templates
    blocker_templates = [
        [(1,5), (12,16), (4,9), (8,13)],
        [(0,4), (6,10), (8,12), (14,18)],
        [(2,6), (4,8), (10,14), (12,16)],
    ]
    for lvl in range(k):
        lo = min(l for l,r in T)
        hi = max(r for l,r in T)
        delta = hi - lo
        pat = copy_patterns[lvl % 4]
        blks = blocker_templates[lvl % 3]
        S = []
        # produce copies according to this level's pattern
        for start in pat:
            off = delta * start - lo
            for (l,r) in T:
                S.append((l+off, r+off))
        # deterministic wave of short intervals to block many active colors
        eps = delta * 0.08
        for w in (3.5, 7.5, 11.5, 15.5):
            st = delta * w - lo
            S.append((st, st + eps))
        # add blockers to couple across copies
        for (a,b) in blks:
            S.append((delta * a, delta * b))
        T = S
    # normalize to small integer grid
    return normalize_to_grid(T)

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()