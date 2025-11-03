# EVOLVE-BLOCK-START

from math import gcd
from functools import reduce

def overlaps(a, b):
    """Open-interval overlap test: return True iff intervals overlap."""
    (l1, r1), (l2, r2) = a, b
    return max(l1, l2) < min(r1, r2)

def firstfit_colors(intervals):
    """
    Simulate FirstFit coloring on the given arrival order.
    Return total number of colors used.
    """
    colors = []  # list of color classes; each is a list of intervals in arrival order
    for iv in intervals:
        placed = False
        for c in colors:
            conflict = False
            # Check conflict within this color class (pairwise disjoint invariant holds,
            # but we conservatively check all to stay robust)
            for u in c:
                if overlaps(u, iv):
                    conflict = True
                    break
            if not conflict:
                c.append(iv)
                placed = True
                break
        if not placed:
            colors.append([iv])
    return len(colors)

def clique_number(intervals):
    """
    Compute omega (maximum number of intervals covering a single point) using sweep.
    For open intervals, endpoints do not contribute to overlap.
    """
    events = []  # (x, type) where type=-1 for right endpoint, +1 for left endpoint
    for (l, r) in intervals:
        if l >= r:
            continue
        events.append((l, +1))
        events.append((r, -1))
    # For open intervals, at the same coordinate handle -1 before +1
    events.sort(key=lambda e: (e[0], 0 if e[1] == -1 else 1))
    cur = 0
    best = 0
    for _, t in events:
        cur += t
        if cur > best:
            best = cur
    return best

def make_copies(T, offsets, delta, lo, center, translation):
    """
    Create 4 translated copies of T according to offsets and translation rule.
    translation in {'left', 'center'}.
    """
    S = []
    for start in offsets:
        if translation == 'left':
            offset = delta * start - lo
        else:  # center-based
            offset = delta * start - center
        for (l, r) in T:
            S.append((l + offset, r + offset))
    return S

def add_blockers(S, blockers, delta, anchor, center):
    """
    Add 4 blockers scaled by delta. Anchor may be 'left' (absolute) or 'center' (center-shifted).
    """
    for (a, b) in blockers:
        if anchor == 'left':
            S.append((delta * a, delta * b))
        else:
            S.append((delta * a - center, delta * b - center))
    return S

def build_pattern(k, base_seed, offsets, blockers, translation, blocker_anchor, schedule='after', permute='id', flip_alt=False):
    """
    Recursively expand the base_seed k times using a multi-copy + blocker scheme with
    engineered arrival orders.

    Args:
      k: recursion depth
      base_seed: initial list of intervals
      offsets: tuple of copy offsets (multipliers of delta)
      blockers: list of 4 pairs (a,b) as multipliers of delta
      translation: 'left' or 'center' anchor for copies
      blocker_anchor: 'left' or 'center' anchor for blockers
      schedule: 'after' | 'before' | 'split' for placing blockers relative to copies
      permute: 'id' | 'rev' | 'zig' reorders the offsets within each level
      flip_alt: if True, alternate reversing the intra-copy arrival order
    """
    def order_offsets(off, mode):
        if mode == 'rev':
            return tuple(reversed(off))
        if mode == 'zig':
            res = []
            i, j = 0, len(off) - 1
            while i <= j:
                res.append(off[i]); i += 1
                if i <= j:
                    res.append(off[j]); j -= 1
            return tuple(res)
        return tuple(off)

    T = list(base_seed)
    for _ in range(k):
        lo = min(l for l, r in T)
        hi = max(r for l, r in T)
        delta = hi - lo
        center = (lo + hi) / 2.0

        # scale blockers
        if blocker_anchor == 'left':
            blk = [(delta * a, delta * b) for (a, b) in blockers]
        else:
            blk = [(delta * a - center, delta * b - center) for (a, b) in blockers]

        # generate copies with engineered order
        offs = order_offsets(offsets, permute)
        copies = []
        for idx, start in enumerate(offs):
            if translation == 'left':
                offset = delta * start - lo
            else:
                offset = delta * start - center
            src = T if not (flip_alt and (idx % 2 == 1)) else list(reversed(T))
            for (l, r) in src:
                copies.append((l + offset, r + offset))

        # Merge according to schedule
        if schedule == 'before':
            S = list(blk) + copies
        elif schedule == 'split':
            mid = len(copies) // 2
            S = copies[:mid] + list(blk) + copies[mid:]
        else:  # 'after'
            S = copies + list(blk)

        T = S
    return T

def normalize_intervals(intervals):
    """
    Normalize to small integer coordinates:
    - scale by 2 to eliminate .5 if produced by center shifts
    - translate so min coordinate is >= 0
    - (optional) divide by global gcd to shrink
    """
    if not intervals:
        return intervals
    # scale by 2 and round (the construction only yields multiples of 0.5)
    scaled = []
    for (l, r) in intervals:
        L = int(round(l * 2))
        R = int(round(r * 2))
        scaled.append((L, R))
    min_coord = min(min(l, r) for l, r in scaled)
    shifted = [(l - min_coord, r - min_coord) for (l, r) in scaled]
    # Keep integers modest—divide by gcd of all endpoints if possible
    vals = []
    for (l, r) in shifted:
        vals.append(abs(l))
        vals.append(abs(r))
    g = 0
    for v in vals:
        g = gcd(g, v)
    if g > 1:
        shrunk = [(l // g, r // g) for (l, r) in shifted]
    else:
        shrunk = shifted
    return shrunk

def evaluate(intervals):
    """
    Compute FF colors and omega, with a small penalty for ridiculously large outputs.
    Return (score, omega, num_colors, n, intervals_normalized)
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
    # soft penalty to prefer smaller n when ratios tie
    score = ratio - 1e-6 * (n / 10000.0)
    return (score, om, cols, n, Tn)

def shrink_optimize_by_ratio(intervals, max_rounds=2):
    """
    Greedy pruning to improve (or preserve) FirstFit/omega ratio.
    Iteratively attempt removing intervals (long-first). Accept a removal if:
      - the ratio strictly increases; or
      - the ratio is unchanged and n decreases (improves the combined score slightly).
    """
    cur = list(intervals)

    def metrics(ints):
        om = clique_number(ints)
        if om <= 0:
            return (0, 0, 0.0)
        cols = firstfit_colors(ints)
        return (cols, om, cols / om)

    cols, om, ratio = metrics(cur)
    if om == 0:
        return cur

    rounds = 0
    improved = True
    while improved and rounds < max_rounds:
        improved = False
        rounds += 1
        # try longer intervals first (tend to be blockers)
        order = sorted(range(len(cur)), key=lambda i: (cur[i][1] - cur[i][0], i), reverse=True)
        for idx in order:
            cand = cur[:idx] + cur[idx + 1:]
            c_cols, c_om, c_ratio = metrics(cand)
            if c_om == 0:
                continue
            if c_ratio > ratio + 1e-12 or (abs(c_ratio - ratio) <= 1e-12 and len(cand) < len(cur)):
                cur = cand
                cols, om, ratio = c_cols, c_om, c_ratio
                improved = True
                break
    return cur

def construct_intervals():
    """
    Build a sequence of open intervals that aims to maximize FirstFit/OPT.
    We sweep several recommended blueprints and pick the best validated candidate.
    Added knobs:
      - schedule of blocker placement within each level
      - permutation of copy offsets (straight/reverse/zigzag)
      - alternating reversal of intra-copy arrival order
      - conservative pruning to improve ratio or reduce size without hurting ratio
    """
    # search space inspired by the provided recommendations
    offsets_set = [
        (2, 6, 10, 14),  # baseline
        (1, 5, 9, 13),
        (3, 7, 11, 15),
        (0, 4, 8, 12),
    ]
    blockers_templates = [
        # Template A: baseline
        ((1, 5), (12, 16), (4, 9), (8, 13)),
        # Template B
        ((0, 4), (6, 10), (8, 12), (14, 18)),
        # Template C
        ((2, 6), (4, 8), (10, 14), (12, 16)),
    ]
    translations = ['left', 'center']  # how copies are positioned
    blocker_anchors = ['left', 'center']  # how blockers are positioned
    depths = [3, 4, 5]  # sweep as recommended
    base_seeds = [
        [(0.0, 1.0)],                      # single seed (classic)
        [(0.0, 1.0), (2.0, 3.0)],          # richer base: two disjoint seeds
    ]
    extra_copies_opts = [0, 1, 2]  # allow up to two extra copies on first level

    schedules = ['after', 'before', 'split']
    permutations = ['id', 'rev', 'zig']
    flip_opts = [False, True]

    best = None  # (score, om, cols, n, intervals_normalized)

    # Enumerate combinations with a guard on worst-case explosion
    for base in base_seeds:
        for k in depths:
            for offsets in offsets_set:
                for extra in extra_copies_opts:
                    # build potentially extended offsets list
                    offs = list(offsets)
                    if extra >= 1:
                        offs.append(offs[-1] + 4)
                    if extra >= 2:
                        offs.append(offs[-1] + 8)
                    offs = tuple(offs)

                    # rough size guard based on branching p=len(offs)
                    p = len(offs)
                    base_sz = len(base)
                    approx = (p ** k) * base_sz + 4 * sum(p ** i for i in range(1, k + 1))
                    if approx > 3000:
                        continue

                    for blockers in blockers_templates:
                        for translation in translations:
                            for anchor in blocker_anchors:
                                for sch in schedules:
                                    for perm in permutations:
                                        for flip in flip_opts:
                                            T = build_pattern(
                                                k=k,
                                                base_seed=base,
                                                offsets=offs,
                                                blockers=blockers,
                                                translation=translation,
                                                blocker_anchor=anchor,
                                                schedule=sch,
                                                permute=perm,
                                                flip_alt=flip,
                                            )
                                            score, om, cols, n, Tn = evaluate(T)
                                            cand = (score, om, cols, n, Tn)
                                            if best is None:
                                                best = cand
                                            else:
                                                # pick better score; tie-break by fewer intervals then larger cols
                                                if cand[0] > best[0] + 1e-9:
                                                    best = cand
                                                elif abs(cand[0] - best[0]) <= 1e-9:
                                                    if cand[3] < best[3]:
                                                        best = cand
                                                    elif cand[3] == best[3] and cand[2] > best[2]:
                                                        best = cand

    # Fallback to baseline if search didn't produce anything
    if best is None:
        # original baseline as a safe default
        k = 4
        T = [(0.0, 1.0)]
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

    # Apply a conservative pruning pass to potentially improve ratio or reduce n
    pruned = shrink_optimize_by_ratio(best[4], max_rounds=2)
    return pruned if pruned else best[4]

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()