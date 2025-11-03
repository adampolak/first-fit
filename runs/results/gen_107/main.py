# EVOLVE-BLOCK-START

def construct_intervals(iterations=4):
    """
    Search small parameter space of recursive expansions, then greedily prune
    to minimize size while preserving FirstFit/opt ≥ target_ratio.
    """
    # simulate FirstFit
    def firstfit_colors(intervals):
        last_end = []
        for (l, r) in intervals:
            placed = False
            for i, e in enumerate(last_end):
                if l >= e:
                    last_end[i] = r
                    placed = True
                    break
            if not placed:
                last_end.append(r)
        return len(last_end)

    # compute clique number via sweep
    def clique_number(intervals):
        events = []
        for l, r in intervals:
            events.append((l, 1))
            events.append((r, -1))
        events.sort(key=lambda x: (x[0], x[1]))
        cur = best = 0
        for _, d in events:
            cur += d
            best = max(best, cur)
        return best

    # build one candidate by recursion
    def make_candidate(depth, starts, extra_first):
        T = [(0.0, 1.0)]
        for level in range(depth):
            lo = min(l for l, r in T)
            hi = max(r for l, r in T)
            delta = hi - lo
            # choose offsets
            offs = list(starts)
            if extra_first and level == 0:
                offs.append(max(starts) + (starts[1]-starts[0]))
            # generate copies
            S = []
            for s in offs:
                off = delta * s - lo
                for (l, r) in T:
                    S.append((l + off, r + off))
            # canonical blockers
            S += [
                (delta*1,  delta*5),
                (delta*12, delta*16),
                (delta*4,  delta*9),
                (delta*8,  delta*13)
            ]
            T = S
        return T

    target_ratio = 2.6
    # parameter grid
    depths = [max(2, iterations-1), iterations]
    start_sets = [
        (2, 6, 10, 14),
        (1, 5, 9, 13),
        (2.5, 6.5, 10.5, 14.5)
    ]
    extras = [False, True]

    best_inst = None
    best_n = float('inf')

    # search
    for d in depths:
        for starts in start_sets:
            for ex in extras:
                inst = make_candidate(d, starts, ex)
                om = clique_number(inst)
                if om == 0:
                    continue
                ff = firstfit_colors(inst)
                if ff/om >= target_ratio:
                    if len(inst) < best_n:
                        best_n = len(inst)
                        best_inst = inst

    # fallback
    if best_inst is None:
        best_inst = make_candidate(iterations, start_sets[0], False)

    # greedy prune by decreasing length
    inst = list(best_inst)
    inst.sort(key=lambda seg: seg[1]-seg[0], reverse=True)
    cur = inst[:]
    improved = True
    while improved:
        improved = False
        for i in range(len(cur)):
            cand = cur[:i] + cur[i+1:]
            om = clique_number(cand)
            if om == 0:
                continue
            ff = firstfit_colors(cand)
            if ff/om >= target_ratio:
                cur = cand
                improved = True
                break
    # normalize to integer grid
    endpoints = sorted({x for seg in cur for x in seg})
    coord = {}
    c = 0
    for x in endpoints:
        coord[x] = c
        c += 2
    normalized = [(coord[l], coord[r]) for (l, r) in cur]
    return normalized

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()