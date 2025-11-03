# EVOLVE-BLOCK-START

import random
from math import gcd
from functools import reduce

# Deterministic randomness for reproducibility
RNG = random.Random(42)

def overlaps(a, b):
    (l1, r1), (l2, r2) = a, b
    return max(l1, l2) < min(r1, r2)

def firstfit_colors(intervals):
    """
    Fast FirstFit coloring using the invariant that intervals assigned to a color
    are non-overlapping, so it's enough to track the right endpoint of the
    last interval assigned to that color (in arrival order).
    """
    last_end = []  # last_end[c-1] = right endpoint of last interval in color c
    for (l, r) in intervals:
        placed = False
        for i, le in enumerate(last_end):
            # if new interval starts after the last interval ended, we can place it
            if l >= le:
                last_end[i] = r
                placed = True
                break
        if not placed:
            last_end.append(r)
    return len(last_end)

def clique_number(intervals):
    """
    Sweep-line for open intervals. For ties, handle right endpoints before left endpoints.
    """
    events = []
    for (l, r) in intervals:
        if l >= r:
            continue
        events.append((l, +1))
        events.append((r, -1))
    events.sort(key=lambda e: (e[0], 0 if e[1] == -1 else 1))
    cur = 0
    best = 0
    for _, t in events:
        cur += t
        if cur > best:
            best = cur
    return best

def normalize_intervals(intervals):
    """
    Convert coordinates to small integers deterministically:
    - multiply by 2 (to remove .5 from center-based constructions)
    - round to nearest integer
    - translate so min coordinate is 0
    - divide by gcd of endpoints to shrink
    """
    if not intervals:
        return []
    scaled = []
    for (l, r) in intervals:
        L = int(round(l * 2))
        R = int(round(r * 2))
        scaled.append((L, R))
    min_coord = min(min(l, r) for l, r in scaled)
    shifted = [(l - min_coord, r - min_coord) for (l, r) in scaled]
    vals = []
    for (l, r) in shifted:
        vals.append(abs(l))
        vals.append(abs(r))
    g = 0
    for v in vals:
        g = gcd(g, v)
    if g > 1:
        return [(l // g, r // g) for (l, r) in shifted]
    else:
        return shifted

def make_copies(T, offsets, delta, lo, center, translation):
    S = []
    for start in offsets:
        if translation == 'left':
            offset = delta * start - lo
        else:
            offset = delta * start - center
        for (l, r) in T:
            S.append((l + offset, r + offset))
    return S

def add_blockers(S, blockers, delta, anchor, center):
    for (a, b) in blockers:
        if anchor == 'left':
            S.append((delta * a, delta * b))
        else:
            S.append((delta * a - center, delta * b - center))
    return S

def build_pattern(k, base_seed, offsets, blockers, translation, blocker_anchor, extra_copies=0):
    """
    Recursively expand base_seed k times using the 4-copy + 4-blocker scheme.
    extra_copies allows adding occasional extra copies on first level (diversify).
    """
    T = list(base_seed)
    for i in range(k):
        lo = min(l for l, r in T)
        hi = max(r for l, r in T)
        delta = hi - lo
        center = (lo + hi) / 2.0

        # choose number of copies: usually len(offsets), optionally +extra at first level
        offs = list(offsets)
        if i == 0 and extra_copies > 0:
            # append a few extra offsets near the end
            for t in range(extra_copies):
                offs.append(offsets[-1] + (t + 1) * 4)

        S = make_copies(T, offs, delta, lo, center, translation)
        S = add_blockers(S, blockers, delta, blocker_anchor, center)
        T = S
    # arrival order: keep as produced (copies then blockers each level)
    return T

def evaluate(intervals):
    """
    Compute normalized intervals, omega, FirstFit colors and score.
    Slight penalty for large n to prefer compact witnesses.
    """
    Tn = normalize_intervals(intervals)
    n = len(Tn)
    if n == 0:
        return (-1.0, 0, 0, n, Tn)
    om = clique_number(Tn)
    if om == 0:
        return (-1.0, 0, 0, n, Tn)
    cols = firstfit_colors(Tn)
    ratio = cols / om
    score = ratio - 1e-6 * (n / 10000.0)
    return (score, om, cols, n, Tn)

def shrink_intervals(intervals, max_trials=2000):
    """
    Try to remove redundant intervals (one by one) while preserving the same
    (cols, om) witness. Greedy randomized order.
    """
    T = list(intervals)
    base_score, base_om, base_cols, _, _ = evaluate(T)
    if len(T) <= 8:
        return T  # already small
    trials = 0
    # candidate indices shuffled
    indices = list(range(len(T)))
    RNG.shuffle(indices)
    i = 0
    while trials < max_trials and i < len(indices):
        idx = indices[i]
        trials += 1
        cand = T[:idx] + T[idx+1:]
        s, om, cols, n, _ = evaluate(cand)
        # accept removal if ratio (cols/om) not decreased and we didn't increase omega
        if om == base_om and cols >= base_cols:
            T = cand
            base_cols = cols
            base_score = s
            # reset search order a bit to attempt further pruning
            indices = list(range(len(T)))
            RNG.shuffle(indices)
            i = 0
            continue
        i += 1
    return T

def random_mutation(params):
    """
    Mutate a parameter dict producing a new one.
    params keys: k, offsets(tuple), blockers(list of pairs), translation, anchor, base_seed_id, extra_copies
    """
    new = dict(params)
    # mutate depth with small prob
    if RNG.random() < 0.2:
        new['k'] = max(2, min(6, new['k'] + RNG.choice([-1, 1])))
    # jitter offsets
    offs = list(new['offsets'])
    for j in range(len(offs)):
        if RNG.random() < 0.5:
            offs[j] = max(0, offs[j] + RNG.choice([-2, -1, 0, 1, 2]))
    # ensure sorted-ish and unique
    offs = sorted(dict.fromkeys(offs))
    # if too few offsets, append some
    while len(offs) < 4:
        offs.append(offs[-1] + 4)
    new['offsets'] = tuple(offs[:4])
    # mutate blockers
    block = [list(b) for b in new['blockers']]
    for j in range(len(block)):
        if RNG.random() < 0.5:
            a, b = block[j]
            # tweak endpoints slightly
            a = max(0, a + RNG.choice([-2, -1, 0, 1, 2]))
            b = max(a+1, b + RNG.choice([-2, -1, 0, 1, 2]))
            block[j] = [a, b]
    # ensure monotone ordering for readability
    block = sorted(block, key=lambda x: (x[0], x[1]))
    new['blockers'] = [tuple(x) for x in block[:4]]
    # toggle translation/anchor sometimes
    if RNG.random() < 0.1:
        new['translation'] = 'center' if new['translation'] == 'left' else 'left'
    if RNG.random() < 0.1:
        new['anchor'] = 'center' if new['anchor'] == 'left' else 'left'
    # change base seed occasionally
    if RNG.random() < 0.1:
        new['base_seed_id'] = 1 - new['base_seed_id']
    # tweak extra_copies
    if RNG.random() < 0.15:
        new['extra_copies'] = max(0, min(3, new['extra_copies'] + RNG.choice([-1, 0, 1])))
    return new

def construct_intervals():
    """
    Stochastic search over the 4-copy+4-blocker blueprint space to find a
    sequence that maximizes FirstFit / omega. Returns normalized intervals.
    """
    # initial templates including the known-good baseline
    initial_offsets = [
        (2, 6, 10, 14),
        (1, 5, 9, 13),
        (3, 7, 11, 15),
        (0, 4, 8, 12),
    ]
    initial_blockers = [
        ((1, 5), (12, 16), (4, 9), (8, 13)),
        ((0, 4), (6, 10), (8, 12), (14, 18)),
        ((2, 6), (4, 8), (10, 14), (12, 16)),
        ((1, 6), (7, 11), (9, 13), (14, 19)),
    ]
    base_seeds_list = [
        [(0.0, 1.0)],
        [(0.0, 1.0), (2.0, 3.0)]
    ]
    # initial parameter
    params = {
        'k': 4,
        'offsets': initial_offsets[0],
        'blockers': list(initial_blockers[0]),
        'translation': 'left',
        'anchor': 'left',
        'base_seed_id': 0,
        'extra_copies': 0
    }

    # evaluate initial pool of seeds
    best_cand = None
    # seed candidate pool with several deterministic variants
    pool = []
    for offs in initial_offsets:
        for blk in initial_blockers:
            for k in (3,4,5):
                pool.append({
                    'k': k,
                    'offsets': offs,
                    'blockers': list(blk),
                    'translation': 'left',
                    'anchor': 'left',
                    'base_seed_id': 0,
                    'extra_copies': 0
                })
    # add a few center-anchor variants
    pool.append({'k':4, 'offsets':(2,6,10,14), 'blockers':list(initial_blockers[0]), 'translation':'center', 'anchor':'center', 'base_seed_id':0, 'extra_copies':0})

    # deterministic search budget
    iterations = 1200
    # limit sizes to keep runtime modest
    max_size = 2200

    # evaluate initial pool
    for p in pool:
        base = base_seeds_list[p['base_seed_id']]
        # approximate max size
        est_size = (len(base) + 4 + p.get('extra_copies',0)) * (len(p['offsets']) ** p['k'] if p['k'] <= 5 else 4**p['k'])
        # filter by a rough cap to avoid explosion
        if est_size > max_size * 10:
            continue
        T = build_pattern(
            k=p['k'],
            base_seed=base,
            offsets=p['offsets'],
            blockers=p['blockers'],
            translation=p['translation'],
            blocker_anchor=p['anchor'],
            extra_copies=p.get('extra_copies', 0)
        )
        score, om, cols, n, Tn = evaluate(T)
        cand = (score, om, cols, n, Tn, p)
        if best_cand is None or cand[0] > best_cand[0]:
            best_cand = cand

    # random local search / hillclimb
    current = {
        'k': 4,
        'offsets': (2, 6, 10, 14),
        'blockers': [(1, 5), (12, 16), (4, 9), (8, 13)],
        'translation': 'left',
        'anchor': 'left',
        'base_seed_id': 0,
        'extra_copies': 0
    }
    # initialize current with best found so far (if exists)
    if best_cand is not None:
        current = best_cand[5]

    current_score = None
    if best_cand is not None:
        current_score = best_cand[0]
    else:
        # evaluate current to have baseline
        T = build_pattern(
            k=current['k'],
            base_seed=base_seeds_list[current['base_seed_id']],
            offsets=current['offsets'],
            blockers=current['blockers'],
            translation=current['translation'],
            blocker_anchor=current['anchor'],
            extra_copies=current.get('extra_copies',0)
        )
        current_score = evaluate(T)[0]

    for it in range(iterations):
        cand_params = random_mutation(current)
        base = base_seeds_list[cand_params['base_seed_id']]
        # quick size estimate: (approx)
        approx_size = (len(base) + 4 + cand_params.get('extra_copies',0)) * (len(cand_params['offsets']) ** min(cand_params['k'],5))
        if approx_size > max_size:
            # skip overly large candidate
            continue
        T = build_pattern(
            k=cand_params['k'],
            base_seed=base,
            offsets=cand_params['offsets'],
            blockers=cand_params['blockers'],
            translation=cand_params['translation'],
            blocker_anchor=cand_params['anchor'],
            extra_copies=cand_params.get('extra_copies', 0)
        )
        s, om, cols, n, Tn = evaluate(T)
        # prefer higher ratio; break ties by fewer intervals then more colors
        better = False
        if best_cand is None or s > best_cand[0] + 1e-9:
            better = True
        elif best_cand is not None and abs(s - best_cand[0]) <= 1e-9:
            if n < best_cand[3]:
                better = True
            elif n == best_cand[3] and cols > best_cand[2]:
                better = True

        if better:
            best_cand = (s, om, cols, n, Tn, cand_params)
            current = cand_params
            current_score = s
        else:
            # accept occasionally for exploration (simulated annealing-like)
            if RNG.random() < 0.01:
                current = cand_params
                current_score = s

    # If nothing found, fall back to baseline
    if best_cand is None:
        # build baseline
        T = [(0.0, 1.0)]
        k = 4
        for _ in range(k):
            lo = min(l for l, r in T)
            hi = max(r for l, r in T)
            delta = hi - lo
            S = []
            for start in (2, 6, 10, 14):
                offset = delta * start - lo
                for (l, r) in T:
                    S.append((l + offset, r + offset))
            S += [
                (delta * 1,  delta * 5),
                (delta * 12, delta * 16),
                (delta * 4,  delta * 9),
                (delta * 8,  delta * 13),
            ]
            T = S
        return normalize_intervals(T)

    # Post-process: try to shrink the witness while preserving score
    best_intervals = best_cand[4]
    # convert back to float construction for shrinker input: expand normalized ints to floats
    # We'll feed shrinker normalized integer list (it uses evaluate which normalizes again) — fine.
    shrunk = shrink_intervals(best_intervals, max_trials=2000)
    return normalize_intervals(shrunk)

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()