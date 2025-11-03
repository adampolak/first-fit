# EVOLVE-BLOCK-START

def construct_intervals(iterations: int = 4):
    """
    Search for the best adversarial interval sequence by varying recursion depth,
    blocker template, and optional extra copy. Return the normalized list of
    integer intervals that maximizes FirstFit/omega ratio.
    """
    # --- Helpers ---
    def normalize(intervals):
        pts = sorted({x for seg in intervals for x in seg})
        m = {x: i*2 for i, x in enumerate(pts)}
        return [(m[l], m[r]) for l, r in intervals]

    def clique_number(intervals):
        events = []
        for l, r in intervals:
            if l < r:
                events.append((l, +1))
                events.append((r, -1))
        events.sort(key=lambda e: (e[0], 0 if e[1] == -1 else 1))
        cur = best = 0
        for _, t in events:
            cur += t
            if cur > best: best = cur
        return best

    def firstfit(intervals):
        # track end of last interval in each color
        ends = []
        for l, r in intervals:
            placed = False
            for i in range(len(ends)):
                if l >= ends[i]:
                    ends[i] = r
                    placed = True
                    break
            if not placed:
                ends.append(r)
        return len(ends)

    def build_pattern(k, blockers, extra_first):
        T = [(0.0, 1.0)]
        base_offsets = (2, 6, 10, 14)
        for i in range(k):
            lo = min(l for l, _ in T)
            hi = max(r for _, r in T)
            delta = hi - lo
            # decide copy offsets
            offs = list(base_offsets)
            if extra_first and i == 0:
                offs.append(18)
            S = []
            # place copies
            for start in offs:
                off = delta*start - lo
                for l, r in T:
                    S.append((l + off, r + off))
            # place blockers for this level
            for a, b in blockers:
                S.append((delta*a, delta*b))
            T = S
        return T

    # Alternative blocker templates (all scaled by delta each level)
    templates = [
        ((1,5), (12,16), (4,9), (8,13)),   # baseline
        ((0,4), (6,10), (8,12), (14,18)),  # template B
        ((1,6), (11,16), (4,8), (10,14)),  # template C
    ]

    best = None  # (ratio, n, cols, omega, normalized_intervals)
    depths = [max(1, iterations-1), iterations, iterations+1]
    for k in depths:
        for tpl in templates:
            for extra in (False, True):
                raw = build_pattern(k, tpl, extra)
                norm = normalize(raw)
                om = clique_number(norm)
                if om == 0:
                    continue
                cols = firstfit(norm)
                ratio = cols / om
                n = len(norm)
                cand = (ratio, n, cols, om, norm)
                if best is None:
                    best = cand
                else:
                    br, bn, bc, bo, _ = best
                    # prefer higher ratio, then fewer intervals, then more colors
                    if ratio > br + 1e-12 or (
                       abs(ratio - br) <= 1e-12 and
                       (n < bn or (n == bn and cols > bc))):
                        best = cand

    # Fallback to classic if search fails
    if best is None:
        T = [(0.0,1.0)]
        for _ in range(iterations):
            lo = min(l for l, _ in T)
            hi = max(r for _, r in T)
            delta = hi - lo
            S = []
            for start in (2,6,10,14):
                off = delta*start - lo
                for l, r in T:
                    S.append((l+off, r+off))
            S += [(delta*1,delta*5),(delta*12,delta*16),(delta*4,delta*9),(delta*8,delta*13)]
            T = S
        return normalize(T)

    # return the best normalized witness
    return best[4]

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()