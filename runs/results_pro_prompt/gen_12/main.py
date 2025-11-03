# EVOLVE-BLOCK-START

import math

def first_fit_color(intervals):
    """Simulate FirstFit on the given interval sequence."""
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

def clique_bound(intervals):
    """Compute the maximum clique size (OPT) via sweep-line on open intervals."""
    events = []
    for l, r in intervals:
        events.append((l, 1))
        events.append((r, -1))
    # process end (-1) before start (+1) at ties
    events.sort(key=lambda e: (e[0], e[1]))
    cur = best = 0
    for _, delta in events:
        cur += delta
        best = max(best, cur)
    return best

def build_fractal(base, depth, starts, caps):
    """Generate a fractal interval set with given start‐shifts and cap intervals."""
    T = list(base)
    for _ in range(depth):
        lo = min(l for l, _ in T)
        hi = max(r for _, r in T)
        delta = hi - lo
        S = []
        # scaled copies
        for s in starts:
            off = delta * s
            S += [(off + l - lo, off + r - lo) for l, r in T]
        # cap intervals
        for a, b in caps:
            S.append((delta * a, delta * b))
        T = S
    return T

def block_hybrid_order(seq):
    """Split seq into 4 blocks and reorder each with a different heuristic."""
    n = len(seq)
    bs = math.ceil(n / 4)
    blocks = [seq[i*bs:(i+1)*bs] for i in range(4)]
    b0 = blocks[0]  # identity
    b1 = list(reversed(blocks[1]))
    b2 = sorted(blocks[2], key=lambda iv: (iv[1]-iv[0], iv[0]))
    b3 = sorted(blocks[3], key=lambda iv: (-(iv[1]-iv[0]), iv[0]))
    return b0 + b1 + b2 + b3

def construct_intervals():
    """
    Build multiple fractal candidates (different anchors + a hybrid),
    test them under several orderings, and return the best FirstFit/OPT.
    """
    base = [(0,1)]
    depth = 4
    # three anchor patterns
    seeds = [
        {'starts': (2,6,10,14), 'caps': [(1,5),(4,9),(8,13),(12,16)]},
        {'starts': (3,7,11,15), 'caps': [(2,6),(5,10),(9,14),(13,17)]},
        {'starts': (4,8,12,16), 'caps': [(3,7),(6,11),(10,15),(14,18)]},
    ]
    candidates = []
    for i, sd in enumerate(seeds):
        seq = build_fractal(base, depth, sd['starts'], sd['caps'])
        candidates.append(('seed{}'.format(i), seq))
    # add one block-hybrid variant of seed0
    candidates.append(('seed0_hybrid', block_hybrid_order(candidates[0][1])))

    # define a suite of ordering heuristics
    orders = [
        ('identity',    lambda x: x),
        ('reversed',    lambda x: list(reversed(x))),
        ('short_first', lambda x: sorted(x, key=lambda iv: (iv[1]-iv[0], iv[0]))),
        ('long_first',  lambda x: sorted(x, key=lambda iv: (-(iv[1]-iv[0]), iv[0]))),
        ('left_first',  lambda x: sorted(x, key=lambda iv: (iv[0], iv[1]))),
        ('right_first', lambda x: sorted(x, key=lambda iv: (-iv[0], -iv[1]))),
    ]

    best_seq = None
    best_ratio = -1.0

    # evaluate every candidate under every ordering
    for cname, seq in candidates:
        for oname, ofunc in orders:
            seq_o = ofunc(seq)
            alg = first_fit_color(seq_o)
            opt = clique_bound(seq_o)
            if opt <= 0:
                continue
            ratio = alg / opt
            if ratio > best_ratio:
                best_ratio = ratio
                best_seq = seq_o

    return best_seq

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()