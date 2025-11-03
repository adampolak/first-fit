# EVOLVE-BLOCK-START

def construct_intervals():
    """
    Search over multiple tiling offsets, blocker templates, and recursion depths.
    For each candidate, build the interval sequence, simulate FirstFit, compute OPT,
    and select the best ratio instance. Returns the normalized integer‐grid intervals
    of that best instance.
    """

    # Definitions of template variants
    start_sets = [
        (2, 6, 10, 14),
        (1, 5, 9, 13),
        (3, 7, 11, 15),
        (0, 4, 8, 12),
    ]
    blocker_templates = {
        "A": [(1,5), (12,16), (4,9), (8,13)],
        "B": [(0,4), (6,10), (8,12), (14,18)],
        "C": [(2,6), (4,8), (10,14), (12,16)],
    }
    best = {
        "ratio": 0.0,
        "intervals": [(0.0,1.0)],
    }

    def build_sequence(starts, blockers, k):
        T = [(0.0,1.0)]
        for _ in range(k):
            lo = min(l for l,r in T)
            hi = max(r for l,r in T)
            delta = hi - lo
            S = []
            # 4 copies
            for s in starts:
                off = delta * s - lo
                for (l,r) in T:
                    S.append((l+off, r+off))
            # blockers, scaled by delta
            for (b0,b1) in blockers:
                S.append((delta*b0, delta*b1))
            T = S
        return T

    def simulate_FF(intervals):
        max_end = []  # max_end[c] = end of last interval assigned color c
        for (l,r) in intervals:
            # try assign to existing color
            assigned = False
            for c in range(len(max_end)):
                if max_end[c] <= l:
                    max_end[c] = r
                    assigned = True
                    break
            if not assigned:
                max_end.append(r)
        return len(max_end)

    def compute_clique(intervals):
        events = []
        for (l,r) in intervals:
            events.append((l, +1))
            events.append((r, -1))
        events.sort(key=lambda x:(x[0], -x[1]))
        cur = mx = 0
        for _,d in events:
            cur += d
            if cur>mx: mx=cur
        return mx

    # Search
    for k in (3,4,5):
        for starts in start_sets:
            for name, blockers in blocker_templates.items():
                seq = build_sequence(starts, blockers, k)
                ff = simulate_FF(seq)
                opt = compute_clique(seq)
                ratio = ff / opt
                # keep best
                if ratio > best["ratio"] + 1e-9:
                    best = {"ratio": ratio, "intervals": seq, "ff":ff, "opt":opt, "desc":(k,starts,name)}

    # Normalize endpoints to integer grid
    T = best["intervals"]
    eps = 1e-9
    pts = sorted(set(x for seg in T for x in seg))
    coord = {}
    cur = 0
    for x in pts:
        coord[x] = cur
        cur += 2
    normalized = [(coord[l], coord[r]) for (l,r) in T]
    # Optionally log chosen template:
    # print("Chosen template:", best["desc"], "FF=",best["ff"],"OPT=",best["opt"],"ratio=",best["ratio"])
    return normalized

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()