# EVOLVE-BLOCK-START

from bisect import bisect_left
from math import gcd

# --------- Basic interval utilities (open intervals) ---------

def _overlap(a, b):
    (l1, r1), (l2, r2) = a, b
    return max(l1, l2) < min(r1, r2)  # open intervals

def _normalize_grid(intervals):
    """
    Map unique endpoints to 0,2,4,... preserving order.
    """
    if not intervals:
        return []
    endpoints = sorted({x for seg in intervals for x in seg})
    coord = {v: i*2 for i, v in enumerate(endpoints)}
    return [(coord[l], coord[r]) for (l, r) in intervals]

def _normalize_grid_gcd(intervals):
    """
    Normalize to 0,2,4,... then divide all coordinates by gcd to shrink.
    """
    T = _normalize_grid(intervals)
    g = 0
    for (l, r) in T:
        g = gcd(g, abs(l))
        g = gcd(g, abs(r))
    if g > 1:
        T = [(l//g, r//g) for (l, r) in T]
    return T

def _omega_open(intervals):
    """
    Maximum number of intervals covering a single point (open intervals).
    Sweep events; process rights before lefts on ties.
    """
    events = []
    for (l, r) in intervals:
        if l < r:
            events.append((l, +1))
            events.append((r, -1))
    events.sort(key=lambda e: (e[0], 0 if e[1] == -1 else 1))
    cur = best = 0
    for _, t in events:
        cur += t
        if cur > best:
            best = cur
    return best

# --------- Accurate FirstFit (per-color non-overlap via neighbors) ---------

def _fits_in_color(iv, color):
    """
    Check if interval iv=(l,r) can be inserted into color without overlap.
    color is {'ivals': [(l,r),...], 'lefts':[l,...]} with ivals sorted by left.
    """
    l, r = iv
    lefts = color['lefts']
    pos = bisect_left(lefts, l)
    # predecessor
    if pos > 0:
        pl, pr = color['ivals'][pos-1]
        if pr > l:  # open intervals => pr > l means overlap
            return False
    # successor
    if pos < len(color['ivals']):
        nl, nr = color['ivals'][pos]
        if nl < r:  # open intervals => nl < r means overlap
            return False
    return True

def _place_in_firstfit(iv, colors):
    """
    Place iv into first color it fits (FirstFit). Returns color index used.
    Mutates colors in place.
    """
    for i, color in enumerate(colors):
        if _fits_in_color(iv, color):
            pos = bisect_left(color['lefts'], iv[0])
            color['ivals'].insert(pos, iv)
            color['lefts'].insert(pos, iv[0])
            return i
    # new color
    colors.append({'ivals':[iv], 'lefts':[iv[0]]})
    return len(colors) - 1

def _firstfit_colors(intervals):
    """
    Return number of colors used by FirstFit on the given ordered list.
    """
    colors = []
    for iv in intervals:
        _place_in_firstfit(iv, colors)
    return len(colors)

# --------- Recursive builders (multi-tiling + blocker templates) ---------

_OFFSETS_SETS = {
    'A': (2, 6, 10, 14),
    'B': (1, 5, 9, 13),
    'C': (3, 7, 11, 15),
    'D': (0, 4, 8, 12),
}

# Three deterministic blocker templates (as multipliers of delta)
_BLOCKERS_TEMPLATES = {
    'A': ((1,5), (12,16), (4,9), (8,13)),                  # canonical
    'B': ((0,4), (7,11), (3,7), (10,14)),                  # template B
    'C': ((1,5), (9,13), (5,9), (12,16)),                  # template C
}

def _expand_once(T, offsets, blockers, schedule='after', translation='left'):
    """
    One expand step: place copies of T at given offsets and append blockers.
    schedule: 'after' | 'before' | 'split'
    translation: 'left' uses offset*delta - lo; 'center' uses offset*delta - center
    """
    lo = min(l for l, r in T)
    hi = max(r for l, r in T)
    delta = hi - lo
    center = (lo + hi) / 2.0

    def make_copies(from_T, offs):
        S = []
        for start in offs:
            if translation == 'left':
                shift = delta * start - lo
            else:
                shift = delta * start - center
            S.extend(((l+shift, r+shift) for (l, r) in from_T))
        return S

    if translation == 'left':
        blk = [(delta*a, delta*b) for (a,b) in blockers]
    else:
        blk = [(delta*a - center, delta*b - center) for (a,b) in blockers]

    if schedule == 'after':
        S = make_copies(T, offsets) + blk
    elif schedule == 'before':
        S = list(blk) + make_copies(T, offsets)
    else:  # split
        h = len(offsets) // 2
        S = make_copies(T, offsets[:h]) + list(blk) + make_copies(T, offsets[h:])
    return S

def _build_blueprint(depth, offset_cycle, blocker_cycle, schedule='after',
                     extra_first=False, translation='left', norm_schedule='final'):
    """
    Build a recursive multi-tiling instance:
    - offset_cycle: string over {'A','B','C','D'} (e.g., 'ABCD' or 'AAAA'),
      selecting offsets per level cyclically.
    - blocker_cycle: string over {'A','B','C'} selecting blockers per level cyclically.
    - schedule: 'after'|'before'|'split' for per-level ordering.
    - extra_first: add a 5th copy offset=18 at level 0 to diversify geometry.
    - translation: 'left'|'center' anchor.
    - norm_schedule: 'final' | 'per_level' | 'two_step' (normalize after every two levels).
    """
    T = [(0.0, 1.0)]
    oc = offset_cycle
    bc = blocker_cycle
    for lvl in range(depth):
        offs = list(_OFFSETS_SETS[oc[lvl % len(oc)]])
        if extra_first and lvl == 0 and 18 not in offs:
            offs.append(18)
        blks = _BLOCKERS_TEMPLATES[bc[lvl % len(bc)]]
        T = _expand_once(T, tuple(offs), blks, schedule=schedule, translation=translation)
        if norm_schedule == 'per_level':
            T = _normalize_grid(T)
        elif norm_schedule == 'two_step' and (lvl % 2 == 1):
            T = _normalize_grid(T)
    return T

# --------- Deterministic adversarial scheduler (novel) ---------

def _midpoint(iv):
    l, r = iv
    return (l + r) / 2.0

def _length(iv):
    return iv[1] - iv[0]

def _covers_point(iv, p):
    l, r = iv
    return l < p and p < r

def _greedy_adversarial_order(intervals, seed_point=None, seed_k=None, pool_k=64):
    """
    Reorder intervals to maximize FirstFit colors:
    - Seed with up to seed_k long intervals covering a common point (spine),
      which forces an initial color stack.
    - Then iteratively pick from a pool of up to pool_k remaining intervals
      the candidate that yields the largest color index if placed next.
    Deterministic: ties are broken by (assigned_color, length, left, right).
    """
    if not intervals:
        return []

    # Normalize a scratch copy for stable comparisons but keep real coords for final return
    S = list(intervals)

    # Choose a seed point near the global center if not given
    if seed_point is None:
        lo = min(l for l, r in S)
        hi = max(r for l, r in S)
        seed_point = (lo + hi) / 2.0

    # Determine seed_k as a safe bound: do not exceed observed omega of the pool
    base_omega = _omega_open(_normalize_grid(S))
    if seed_k is None:
        seed_k = max(3, min(7, base_omega))  # conservative; typically equals omega

    # Spine: take the longest intervals covering seed_point
    cover = [iv for iv in S if _covers_point(iv, seed_point)]
    cover.sort(key=lambda iv: (-_length(iv), iv[0], iv[1]))
    spine = cover[:seed_k]
    used = set(spine)
    ordered = []
    colors = []

    # Place spine intervals first (they will stack into distinct colors)
    for iv in spine:
        _place_in_firstfit(iv, colors)
        ordered.append(iv)

    # Remaining intervals
    remaining = [iv for iv in S if iv not in used]

    # Greedy loop
    while remaining:
        # Candidate pool: pick up to pool_k longest remaining intervals
        remaining.sort(key=lambda iv: (-_length(iv), iv[0], iv[1]))
        pool = remaining[:pool_k]

        best = None  # (assigned_color, length, left, right, iv, idx)
        best_idx = None
        # Evaluate candidates on current colors
        for idx, iv in enumerate(pool):
            # simulate placement
            assigned = None
            for i, color in enumerate(colors):
                if _fits_in_color(iv, color):
                    assigned = i
                    break
            if assigned is None:
                assigned = len(colors)  # would open a new color

            key = (assigned, _length(iv), -iv[0], -iv[1])  # prefer larger color, longer, right-shifted
            if (best is None) or (key > best[0]):
                best = (key, iv)
                best_idx = idx

        # Fallback if pool evaluated empty (should not happen)
        if best is None:
            iv = remaining.pop(0)
        else:
            iv = pool[best_idx]
            remaining.remove(iv)

        _place_in_firstfit(iv, colors)
        ordered.append(iv)

    return ordered

# --------- Frontier-preserving pruning (deterministic) ---------

def _evaluate(intervals):
    Tn = _normalize_grid(intervals)
    omega = _omega_open(Tn)
    cols = _firstfit_colors(Tn)
    ratio = cols / omega if omega > 0 else 0.0
    return ratio, cols, omega, len(Tn), Tn

def _prune_preserve_frontier(intervals, target_ratio=None):
    """
    Remove intervals if colors/omega does not drop below target_ratio.
    Earliest-added-first removal, one pass then selective repeats.
    """
    cur = list(intervals)
    if not cur:
        return cur

    base_ratio, base_cols, base_omega, _, _ = _evaluate(cur)
    if target_ratio is None:
        target_ratio = base_ratio

    changed = True
    while changed:
        changed = False
        for i in range(len(cur)):
            cand = cur[:i] + cur[i+1:]
            r, c, o, _, _ = _evaluate(cand)
            if o == 0:
                continue
            if r >= target_ratio - 1e-12:
                cur = cand
                changed = True
                break
    return cur

# --------- High-level construction and selection ---------

def construct_intervals(iterations=4, normalize=True):
    """
    Deterministic exploration + adversarial scheduling + pruning.
    Returns a list of (l,r) open intervals in arrival order.
    """

    # Candidate parameter grids (compact, diverse)
    depths = sorted(set([max(2, iterations - 1), max(2, iterations), min(5, max(2, iterations + 1))]))
    offset_cycles = [
        'AAAA', 'BBBB', 'CCCC', 'DDDD',
        'ABCD', 'BCDA', 'CDAB', 'DABC',
        'ACAC', 'BDBD'
    ]
    blocker_cycles = [
        'A', 'B', 'C', 'AB', 'AC', 'BC', 'ABC'
    ]
    schedules = ['after', 'before', 'split']
    translations = ['left', 'center']
    norm_schedules = ['final', 'per_level', 'two_step']
    extra_first_opts = [False, True]

    best_data = None  # (ratio, cols, omega, n, raw_intervals)
    best_T_raw = None

    # Stage 1: build blueprints and evaluate raw
    for d in depths:
        # limit explosion guard (approximate growth ~ 4^d)
        if d > 5:
            continue
        for oc in offset_cycles:
            for bc in blocker_cycles:
                for sch in schedules:
                    for tr in translations:
                        for ns in norm_schedules:
                            for ef in extra_first_opts:
                                T = _build_blueprint(d, oc, bc, schedule=sch,
                                                     extra_first=ef, translation=tr, norm_schedule=ns)
                                ratio, cols, omega, n, Tn = _evaluate(T)
                                # Prefer higher ratio, then fewer intervals, then higher cols
                                if omega == 0:
                                    continue
                                cand = (ratio, -n, cols, T)
                                if (best_data is None) or (cand > (best_data[0], -best_data[3], best_data[1], best_T_raw)):
                                    best_data = (ratio, cols, omega, n, Tn)
                                    best_T_raw = T

    # Fallback to canonical if search failed
    if best_T_raw is None:
        T = [(0.0, 1.0)]
        for _ in range(max(2, iterations)):
            lo = min(l for l, r in T)
            hi = max(r for l, r in T)
            delta = hi - lo
            S = []
            for start in (2, 6, 10, 14):
                S.extend(((delta*start + l - lo, delta*start + r - lo) for (l, r) in T))
            S += [(delta*1, delta*5), (delta*12, delta*16), (delta*4, delta*9), (delta*8, delta*13)]
            T = S
        best_T_raw = T

    # Stage 2: adversarial scheduling (novel reordering)
    # Determine seed point as center of span; seed_k bounded by observed omega
    lo = min(l for l, r in best_T_raw)
    hi = max(r for l, r in best_T_raw)
    center_point = (lo + hi) / 2.0
    raw_ratio, raw_cols, raw_omega, raw_n, _ = _evaluate(best_T_raw)

    reordered = _greedy_adversarial_order(best_T_raw, seed_point=center_point, seed_k=None, pool_k=64)
    re_ratio, re_cols, re_omega, re_n, _ = _evaluate(reordered)

    # Choose better between raw and reordered
    if (re_ratio > raw_ratio + 1e-12) or (abs(re_ratio - raw_ratio) <= 1e-12 and (re_n < raw_n or re_cols > raw_cols)):
        chosen = reordered
        base_ratio = re_ratio
    else:
        chosen = best_T_raw
        base_ratio = raw_ratio

    # Stage 3: frontier-preserving pruning to reduce footprint without hurting frontier
    pruned = _prune_preserve_frontier(chosen, target_ratio=base_ratio)
    fin_ratio, fin_cols, fin_omega, fin_n, Tn = _evaluate(pruned)

    # Normalization for output
    if normalize:
        return _normalize_grid_gcd(pruned if fin_ratio >= base_ratio - 1e-12 else chosen)
    else:
        return pruned if fin_ratio >= base_ratio - 1e-12 else chosen

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()