# EVOLVE-BLOCK-START
def construct_intervals():
    import random, math
    random.seed(0)

    # FirstFit simulator
    def firstfit_colors(seq):
        end_times = []
        for l, r in seq:
            placed = False
            for i, et in enumerate(end_times):
                if et <= l:
                    end_times[i] = r
                    placed = True
                    break
            if not placed:
                end_times.append(r)
        return len(end_times)

    # Offline clique number
    def clique_number(seq):
        events = []
        for l, r in seq:
            events.append((l, 1))
            events.append((r, -1))
        events.sort(key=lambda e: (e[0], e[1]))
        cur = best = 0
        for _, d in events:
            cur += d
            if cur > best:
                best = cur
        return best

    # Block‐hybrid arrival order
    def block_hybrid(seq):
        n = len(seq); bs = (n + 3) // 4
        blocks = [seq[i*bs:(i+1)*bs] for i in range(4)]
        while len(blocks) < 4:
            blocks.append([])
        return (
            blocks[0]
            + list(reversed(blocks[1]))
            + sorted(blocks[2], key=lambda iv: (iv[1] - iv[0], iv[0]))
            + sorted(blocks[3], key=lambda iv: (-(iv[1] - iv[0]), iv[0]))
        )

    # Greedy: pick interval that will be placed in highest color
    def greedy_assigned(seq, sample, seed):
        rng = random.Random(seed)
        rem = list(seq)
        out = []
        ends = []
        while rem:
            k = min(sample, len(rem))
            idxs = list(range(len(rem))) if k == len(rem) else rng.sample(range(len(rem)), k)
            best_idx, best_col = idxs[0], -1
            for idx in idxs:
                l, r = rem[idx]
                col = next((i for i, et in enumerate(ends) if et <= l), len(ends))
                if col > best_col:
                    best_idx, best_col = idx, col
            l, r = rem.pop(best_idx)
            placed = False
            for i, et in enumerate(ends):
                if et <= l:
                    ends[i] = r
                    placed = True
                    break
            if not placed:
                ends.append(r)
            out.append((l, r))
        return out

    # Greedy: pick interval overlapping most active colors
    def greedy_overlap(seq, sample, seed):
        rng = random.Random(seed)
        rem = list(seq)
        out = []
        ends = []
        while rem:
            k = min(sample, len(rem))
            idxs = list(range(len(rem))) if k == len(rem) else rng.sample(range(len(rem)), k)
            best_idx, best_ov = idxs[0], -1
            for idx in idxs:
                l, r = rem[idx]
                ov = sum(1 for et in ends if et > l)
                if ov > best_ov:
                    best_idx, best_ov = idx, ov
            l, r = rem.pop(best_idx)
            placed = False
            for i, et in enumerate(ends):
                if et <= l:
                    ends[i] = r
                    placed = True
                    break
            if not placed:
                ends.append(r)
            out.append((l, r))
        return out

    # Fractal builder with optional interleaving & caps‐first
    def build_fractal(base, depth, starts, caps, interleave, caps_before):
        T = list(base)
        for _ in range(depth):
            lo = min(l for l, _ in T)
            hi = max(r for _, r in T)
            delta = hi - lo
            S = []
            if caps_before:
                for a, b in caps:
                    S.append((delta * a, delta * b))
            if interleave:
                for (li, ri) in T:
                    for s in starts:
                        off = delta * s
                        S.append((off + li - lo, off + ri - lo))
            else:
                for s in starts:
                    S += [((delta * s + l - lo), (delta * s + r - lo)) for l, r in T]
            if not caps_before:
                for a, b in caps:
                    S.append((delta * a, delta * b))
            T = S
            if len(T) > 2000:
                break
        return [(int(l), int(r)) for l, r in T]

    # Three classic seeds
    seed_patterns = [
        ((2, 6, 10, 14), [(1, 5), (4, 9), (8, 13), (12, 16)]),
        ((3, 7, 11, 15), [(2, 6), (5, 10), (9, 14), (13, 17)]),
        ((4, 8, 12, 16), [(3, 7), (6, 11), (10, 15), (14, 18)]),
    ]

    base = [(0.0, 1.0)]
    best_seq = None
    best_score = -1.0

    # Parameter sweep
    for starts, caps in seed_patterns:
        for depth in (3, 4):
            for inter in (False, True):
                for capb in (False, True):
                    seq = build_fractal(base, depth, starts, caps, inter, capb)
                    if not (1 <= len(seq) <= 2000):
                        continue
                    ω = clique_number(seq)
                    if ω < 1:
                        continue
                    # arrival‐order candidates
                    orders = [
                        seq,
                        list(reversed(seq)),
                        sorted(seq, key=lambda iv: (iv[0], iv[1])),
                        sorted(seq, key=lambda iv: (iv[1] - iv[0], iv[0])),
                        block_hybrid(seq),
                        greedy_assigned(seq, 40, 0),
                        greedy_overlap(seq, 40, 1),
                    ]
                    ratios = [firstfit_colors(o) / ω for o in orders]
                    peak = max(ratios)
                    mean = sum(ratios) / len(ratios)
                    score = peak + 0.1 * mean
                    if score > best_score:
                        best_score = score
                        best_seq = orders[ratios.index(peak)]

    return best_seq or [(0, 1)]
# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()