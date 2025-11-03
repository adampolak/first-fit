# EVOLVE-BLOCK-START

def construct_intervals():
    import math, random
    random.seed(42)

    # FirstFit simulator
    def firstfit_color_count(intervals):
        end_times = []
        for l, r in intervals:
            placed = False
            for i, et in enumerate(end_times):
                if et <= l:
                    end_times[i] = r
                    placed = True
                    break
            if not placed:
                end_times.append(r)
        return len(end_times)

    # Offline optimum via clique number
    def clique_bound(intervals):
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

    # Fractal construction
    def build_fractal(base, depth, starts, caps):
        T = list(base)
        for _ in range(depth):
            lo = min(l for l, _ in T)
            hi = max(r for _, r in T)
            delta = hi - lo
            S = []
            # scaled copies
            for s in starts:
                off = delta * s
                S.extend([(off + l - lo, off + r - lo) for l, r in T])
            # cap intervals
            for a, b in caps:
                S.append((delta * a, delta * b))
            T = S
        return T

    # Block‐hybrid reordering
    def block_hybrid_order(seq):
        n = len(seq)
        if n < 4:
            return list(seq)
        k = 4
        b = n // k
        parts = [seq[i*b:(i+1)*b] for i in range(k-1)]
        parts.append(seq[(k-1)*b:])
        out = []
        out += parts[0]
        if len(parts) > 1:
            out += sorted(parts[1], key=lambda iv: (iv[0], iv[1]))
        if len(parts) > 2:
            out += list(reversed(parts[2]))
        if len(parts) > 3:
            out += sorted(parts[3], key=lambda iv: (iv[1]-iv[0], iv[0]))
        return out

    # Greedy that picks the interval forcing the largest FirstFit color
    def greedy_max_assigned_color(seq, sample_size=80, seed=0):
        rng = random.Random(seed)
        rem = list(seq)
        order = []
        end_times = []
        while rem:
            k = min(sample_size, len(rem))
            idxs = list(range(len(rem))) if k == len(rem) else rng.sample(range(len(rem)), k)
            best_idx = idxs[0]
            best_color = -1
            for idx in idxs:
                l, r = rem[idx]
                assigned = None
                for i, et in enumerate(end_times):
                    if et <= l:
                        assigned = i
                        break
                if assigned is None:
                    assigned = len(end_times)
                if assigned > best_color:
                    best_color = assigned
                    best_idx = idx
            l, r = rem.pop(best_idx)
            placed = False
            for i, et in enumerate(end_times):
                if et <= l:
                    end_times[i] = r
                    placed = True
                    break
            if not placed:
                end_times.append(r)
            order.append((l, r))
        return order

    # Greedy that picks the interval overlapping the most active colors
    def greedy_max_overlap(seq, sample_size=80, seed=1):
        rng = random.Random(seed)
        rem = list(seq)
        order = []
        end_times = []
        while rem:
            k = min(sample_size, len(rem))
            idxs = list(range(len(rem))) if k == len(rem) else rng.sample(range(len(rem)), k)
            best_idx = idxs[0]
            best_ov = -1
            for idx in idxs:
                l, r = rem[idx]
                ov = sum(1 for et in end_times if et > l)
                if ov > best_ov:
                    best_ov = ov
                    best_idx = idx
            l, r = rem.pop(best_idx)
            placed = False
            for i, et in enumerate(end_times):
                if et <= l:
                    end_times[i] = r
                    placed = True
                    break
            if not placed:
                end_times.append(r)
            order.append((l, r))
        return order

    # Seed patterns (starts + caps)
    seeds = [
        {'starts': (2, 6, 10, 14), 'caps': [(1,5),(4,9),(8,13),(12,16)]},
        {'starts': (3, 7, 11, 15), 'caps': [(2,6),(5,10),(9,14),(13,17)]},
        {'starts': (4, 8, 12, 16), 'caps': [(3,7),(6,11),(10,15),(14,18)]},
    ]

    base = [(0.0, 1.0)]
    depth_options = [3, 4]
    max_intervals = 2000

    best_seq = None
    best_ratio = -1.0

    # Explore each seed & depth
    for sd in seeds:
        for depth in depth_options:
            seq_f = build_fractal(base, depth, sd['starts'], sd['caps'])
            if len(seq_f) == 0 or len(seq_f) > max_intervals:
                continue
            # integer endpoints
            seq = [(int(l), int(r)) for l, r in seq_f]
            omega = clique_bound(seq)
            if omega <= 0:
                continue
            # Arrival‐order heuristics
            orderings = [
                ('id',         lambda x: x),
                ('rev',        lambda x: list(reversed(x))),
                ('left',       lambda x: sorted(x, key=lambda iv:(iv[0],iv[1]))),
                ('right',      lambda x: sorted(x, key=lambda iv:(-iv[0],-iv[1]))),
                ('short',      lambda x: sorted(x, key=lambda iv:(iv[1]-iv[0],iv[0]))),
                ('long',       lambda x: sorted(x, key=lambda iv:(-(iv[1]-iv[0]),iv[0]))),
                ('block_hyb',  block_hybrid_order),
                ('g_assn_40',  lambda x: greedy_max_assigned_color(x,40,0)),
                ('g_assn_80',  lambda x: greedy_max_assigned_color(x,80,0)),
                ('g_ovlp_40',  lambda x: greedy_max_overlap(x,40,1)),
                ('g_ovlp_80',  lambda x: greedy_max_overlap(x,80,1)),
            ]
            for name, ofunc in orderings:
                try:
                    seq_o = ofunc(seq)
                except:
                    continue
                alg = firstfit_color_count(seq_o)
                ratio = alg / omega
                if ratio > best_ratio or (abs(ratio-best_ratio)<1e-9 and (best_seq is None or len(seq_o)<len(best_seq))):
                    best_ratio = ratio
                    best_seq = seq_o

    # Fallback to a trivial interval if nothing found
    if best_seq is None:
        best_seq = [(0,1)]
    return best_seq

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()