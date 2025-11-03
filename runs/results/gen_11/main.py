# EVOLVE-BLOCK-START

def construct_intervals():
    """
    Search over a few related variants of the 4-copy / 4-blocker construction
    and return the candidate that maximizes the empirical FirstFit / OPT ratio.

    The search varies:
      - tiling offsets (several symmetric choices)
      - blocker templates (a couple of plausible alternatives)
      - recursion depth (3,4,5)
      - small base seeds (one or two tiny seeds)

    For each candidate we:
      - construct intervals deterministically,
      - simulate FirstFit on that arrival order,
      - compute the offline clique number via sweep-line,
      - measure the ratio and pick the best candidate.

    A hard cap prevents explosion of interval counts.
    """
    # helper: simulate FirstFit for intervals in given order
    def simulate_firstfit(intervals):
        # maintain end time of the last interval assigned to each color
        last_end = []  # last_end[c] = end point of last interval on color c
        for (l, r) in intervals:
            # find smallest color with last_end <= l (non-overlapping)
            assigned = False
            for i in range(len(last_end)):
                if last_end[i] <= l:
                    last_end[i] = r
                    assigned = True
                    break
            if not assigned:
                last_end.append(r)
        return len(last_end)

    # helper: compute clique number (maximum number of intervals covering a point)
    def compute_clique(intervals):
        # events: (x, type) where type = -1 for end, +1 for start.
        # When equal coordinate, process end before start to respect open intervals.
        events = []
        for (l, r) in intervals:
            events.append((l, +1))
            events.append((r, -1))
        events.sort(key=lambda x: (x[0], x[1]))  # ends (-1) come before starts (+1)
        cur = 0
        best = 0
        for pos, typ in events:
            if typ == -1:
                cur -= 1
            else:
                cur += 1
                if cur > best:
                    best = cur
        return best

    # builder: given a seed T, starts tuple and blocker template, build with depth k
    def build_from(seed, starts, blockers, depth, max_intervals=4000):
        T = list(seed)
        for _ in range(depth):
            lo = min(l for l, r in T)
            hi = max(r for l, r in T)
            delta = hi - lo if hi > lo else 1.0
            S = []
            # tiled copies
            for start in starts:
                offset = delta * start - lo
                for (l, r) in T:
                    S.append((l + offset, r + offset))
            # add blockers (each blocker specified as a pair of multiplicative factors)
            for (a, b) in blockers:
                S.append((delta * a, delta * b))
            T = S
            if len(T) > max_intervals:
                # abort this build (too large)
                return None
        return T

    # candidate parameters
    starts_options = [
        (2, 6, 10, 14),
        (1, 5, 9, 13),
        (3, 7, 11, 15),
        (0, 4, 8, 12),
    ]
    blocker_templates = [
        # original (multiplicative factors)
        [(1, 5), (12, 16), (4, 9), (8, 13)],
        # variant B
        [(0, 4), (6, 10), (8, 12), (14, 18)],
        # variant C
        [(2, 6), (4, 8), (10, 14), (12, 16)],
    ]
    seeds = [
        [(0.0, 1.0)],  # single small interval (classic)
        [(0.0, 1.0), (2.0, 3.0)],  # slightly richer base
    ]
    depths = [3, 4, 5]

    best_score = -1.0
    best_intervals = seeds[0]

    # search loop
    for seed in seeds:
        for starts in starts_options:
            for blockers in blocker_templates:
                for depth in depths:
                    intervals = build_from(seed, starts, blockers, depth, max_intervals=4000)
                    if not intervals:
                        continue
                    # simulate
                    ff = simulate_firstfit(intervals)
                    opt = compute_clique(intervals)
                    if opt == 0:
                        continue
                    score = ff / float(opt)
                    # prefer higher score; tie-breaker: fewer intervals then larger ff
                    if score > best_score or (abs(score - best_score) < 1e-9 and (len(intervals) < len(best_intervals) or ff > simulate_firstfit(best_intervals))):
                        best_score = score
                        best_intervals = intervals

    # fallback: if nothing selected (shouldn't happen), return original simple construction
    if not best_intervals:
        return [(0.0, 1.0)]

    return best_intervals

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()