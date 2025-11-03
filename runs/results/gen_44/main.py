# EVOLVE-BLOCK-START

from math import gcd

def overlaps(a, b):
    (l1, r1), (l2, r2) = a, b
    return max(l1, l2) < min(r1, r2)

def firstfit_colors(intervals):
    # efficient FirstFit using last endpoint per color (intervals are open)
    last_end = []
    for (l, r) in intervals:
        placed = False
        for i in range(len(last_end)):
            if l >= last_end[i]:
                last_end[i] = r
                placed = True
                break
        if not placed:
            last_end.append(r)
    return len(last_end)

def clique_number(intervals):
    # sweep-line for open intervals: treat right endpoint before left at ties
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
    # scale by 2 to clear halves, round, shift to nonnegative, divide by gcd
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
    T = list(base_seed)
    for level in range(k):
        lo = min(l for l, r in T)
        hi = max(r for l, r in T)
        delta = hi - lo
        center = (lo + hi) / 2.0
        offs = list(offsets)
        if level == 0 and extra_copies > 0:
            last = offs[-1] if offs else 0
            for t in range(extra_copies):
                offs.append(last + 4 * (t + 1))
        S = make_copies(T, offs, delta, lo, center, translation)
        S = add_blockers(S, blockers, delta, blocker_anchor, center)
        T = S
    return T

def evaluate(intervals):
    Tn = normalize_intervals(intervals)
    n = len(Tn)
    if n == 0:
        return (-1.0, 0, 0, n, Tn)
    om = clique_number(Tn)
    if om == 0:
        return (-1.0, 0, 0, n, Tn)
    cols = firstfit_colors(Tn)
    ratio = cols / om
    # tiny penalty for larger size to prefer compact witnesses
    score = ratio - 1e-6 * (n / 10000.0)
    return (score, om, cols, n, Tn)

def shrink_intervals(intervals, max_trials=5000):
    # greedy deterministic removal from end to start, accepting removals that preserve (om, cols)
    T = list(intervals)
    best_score, best_om, best_cols, _, _ = evaluate(T)
    if len(T) <= 8:
        return T
    trials = 0
    i = 0
    while trials < max_trials and i < len(T):
        idx = len(T) - 1 - i
        trials += 1
        cand = T[:idx] + T[idx+1:]
        s, om, cols, n, _ = evaluate(cand)
        # accept removal if omega unchanged and cols >= previous (don't reduce FirstFit's colors)
        if om == best_om and cols >= best_cols:
            T = cand
            best_cols = cols
            best_score = s
            i = 0
            continue
        i += 1
    return T

def construct_intervals():
    """
    Deterministic search around the 4-copy/4-blocker blueprint.
    Returns normalized list of (l,r) intervals in arrival order.
    """
    # parameter space (kept small & deterministic)
    offsets_list = [
        (2, 6, 10, 14),  # baseline
        (1, 5, 9, 13),
        (3, 7, 11, 15),
        (0, 4, 8, 12),
    ]
    blockers_templates = [
        # baseline and a few variants
        ((1, 5), (12, 16), (4, 9), (8, 13)),
        ((0, 4), (6, 10), (8, 12), (14, 18)),
        ((1, 6), (7, 11), (9, 13), (14, 19)),
    ]
    translations = ['left', 'center']
    anchors = ['left', 'center']
    depths = [3, 4, 5]  # try around proven sweet spot
    base_seeds = [
        [(0.0, 1.0)],
        [(0.0, 1.0), (2.0, 3.0)],
    ]
    extra_opts = [0, 1]  # at most one extra copy on first level
    MAX_INTERVALS = 2200

    best = None
    for base in base_seeds:
        for k in depths:
            for offsets in offsets_list:
                for blockers in blockers_templates:
                    for translation in translations:
                        for anchor in anchors:
                            for extra in extra_opts:
                                # estimate size to avoid explosion
                                copies = len(offsets) + (extra if extra else 0)
                                est = (copies ** k) * len(base) + 4 * (copies ** (k - 1)) * k
                                if est > MAX_INTERVALS:
                                    continue
                                T = build_pattern(
                                    k=k,
                                    base_seed=base,
                                    offsets=offsets,
                                    blockers=blockers,
                                    translation=translation,
                                    blocker_anchor=anchor,
                                    extra_copies=extra
                                )
                                cand_eval = evaluate(T)
                                s, om, cols, n, Tn = cand_eval
                                cand = (s, om, cols, n, Tn, T)
                                if best is None:
                                    best = cand
                                else:
                                    if cand[0] > best[0] + 1e-9:
                                        best = cand
                                    elif abs(cand[0] - best[0]) <= 1e-9:
                                        if cand[3] < best[3]:
                                            best = cand
                                        elif cand[3] == best[3] and cand[2] > best[2]:
                                            best = cand

    # fallback to a modest default if nothing found
    if best is None:
        T = [(0.0, 1.0)]
        k = 4
        for _ in range(k):
            lo = min(l for l, r in T)
            hi = max(r for l, r in T)
            delta = hi - lo
            S = []
            for start in (2, 6, 10, 14):
                S += [(delta * start + l - lo, delta * start + r - lo) for l, r in T]
            S += [
                (delta * 1, delta * 5),
                (delta * 12, delta * 16),
                (delta * 4, delta * 9),
                (delta * 8, delta * 13)
            ]
            T = S
        return normalize_intervals(T)

    # try to shrink the best raw construction
    raw = best[5]
    shrunk = shrink_intervals(raw, max_trials=5000)
    final = normalize_intervals(shrunk)
    # safety: fall back to normalized best if shrink emptied (shouldn't happen)
    if not final:
        return best[4]
    return final

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()