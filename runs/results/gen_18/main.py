# EVOLVE-BLOCK-START

from math import gcd
from functools import reduce

def overlaps(a, b):
    """Open-interval overlap test."""
    (l1, r1), (l2, r2) = a, b
    return max(l1, l2) < min(r1, r2)

def firstfit_colors(intervals):
    """
    Efficient FirstFit simulation: track right endpoint of last interval in each color.
    Since intervals assigned to a color are non-overlapping, it's enough to check start >= last_end.
    """
    last_end = []
    for (l, r) in intervals:
        placed = False
        # Traverse colors in order (FirstFit)
        for i in range(len(last_end)):
            if l >= last_end[i]:
                last_end[i] = r
                placed = True
                break
        if not placed:
            last_end.append(r)
    return len(last_end)

def clique_number(intervals):
    """
    Sweep-line for open intervals. For equal coordinates, handle right endpoints before left endpoints.
    """
    events = []
    for (l, r) in intervals:
        if l >= r:
            continue
        events.append((l, +1))
        events.append((r, -1))
    # sort: by position; for ties, -1 (right) before +1 (left) so open intervals don't touch
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
    Normalize coordinates to small integers:
    - multiply coordinates by 2 (to clear possible .5 centers)
    - round to ints
    - shift so min coord is 0
    - divide by gcd of all endpoints to shrink
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
    """Create translated copies of T according to offsets and translation mode."""
    S = []
    for start in offsets:
        if translation == 'left':
            offset = delta * start - lo
        else:  # 'center'
            offset = delta * start - center
        for (l, r) in T:
            S.append((l + offset, r + offset))
    return S

def add_blockers(S, blockers, delta, anchor, center):
    """Add blockers scaled by delta, anchored either by 'left' or 'center'."""
    for (a, b) in blockers:
        if anchor == 'left':
            S.append((delta * a, delta * b))
        else:
            S.append((delta * a - center, delta * b - center))
    return S

def build_pattern(k, base_seed, offsets, blockers, translation, blocker_anchor, extra_copies=0):
    """
    Expand base_seed k times using the copy + blocker scheme.
    extra_copies (small integer) appends a few additional translated copies on the first level
    to diversify geometry; kept tiny to avoid blowup.
    """
    T = list(base_seed)
    for level in range(k):
        lo = min(l for l, r in T)
        hi = max(r for l, r in T)
        delta = hi - lo
        center = (lo + hi) / 2.0

        offs = list(offsets)
        if level == 0 and extra_copies > 0:
            # Add extra offsets after the baseline ones spaced by 4
            last = offs[-1] if offs else 0
            for t in range(extra_copies):
                offs.append(last + 4 * (t + 1))

        S = make_copies(T, offs, delta, lo, center, translation)
        S = add_blockers(S, blockers, delta, blocker_anchor, center)
        T = S
    return T

def evaluate(intervals):
    """Normalize, compute omega and FirstFit colors, and return score tuple."""
    Tn = normalize_intervals(intervals)
    n = len(Tn)
    if n == 0:
        return (-1.0, 0, 0, n, Tn)
    om = clique_number(Tn)
    if om == 0:
        return (-1.0, 0, 0, n, Tn)
    cols = firstfit_colors(Tn)
    ratio = cols / om
    # small penalty for larger n to prefer compact witnesses
    score = ratio - 1e-6 * (n / 10000.0)
    return (score, om, cols, n, Tn)

def shrink_intervals(intervals, max_trials=5000):
    """
    Greedy attempt to remove intervals while preserving (cols, om).
    Deterministic traversal left-to-right removing if witness remains.
    """
    T = list(intervals)
    base_score, base_om, base_cols, _, _ = evaluate(T)
    if len(T) <= 8:
        return T
    i = 0
    trials = 0
    # deterministic order: try from end to start to remove later redundant blockers first
    while trials < max_trials and i < len(T):
        idx = len(T) - 1 - i
        trials += 1
        cand = T[:idx] + T[idx+1:]
        s, om, cols, n, _ = evaluate(cand)
        # Accept removal if omega unchanged and cols >= base_cols
        if om == base_om and cols >= base_cols:
            T = cand
            base_cols = cols
            base_score = s
            # reset i to try again from end
            i = 0
            continue
        i += 1
    return T

def construct_intervals():
    """
    Deterministic exploration of variants of the 4-copy + 4-blocker blueprint.
    Returns a normalized integer list of intervals (list of (l, r) pairs).
    """
    # Blueprint options inspired by recommendations
    offsets_set = [
        (2, 6, 10, 14),
        (1, 5, 9, 13),
        (3, 7, 11, 15),
        (0, 4, 8, 12),
    ]
    blockers_templates = [
        ((1, 5), (12, 16), (4, 9), (8, 13)),   # baseline A
        ((0, 4), (6, 10), (8, 12), (14, 18)),  # B
        ((2, 6), (4, 8), (10, 14), (12, 16)),  # C
        ((1, 6), (7, 11), (9, 13), (14, 19)),  # variant D
    ]
    translations = ['left', 'center']
    blocker_anchors = ['left', 'center']
    depths = [3, 4, 5]
    base_seeds = [
        [(0.0, 1.0)],                      # single short seed
        [(0.0, 1.0), (2.0, 3.0)],          # richer base
    ]
    extra_copies_opts = [0, 1]  # small diversification

    # size cap to avoid explosion
    MAX_INTERVALS = 2200

    best = None  # best candidate tuple (score, om, cols, n, normalized_intervals, raw_intervals)
    # Enumerate deterministic combinations
    for base_id, base in enumerate(base_seeds):
        base_size = len(base)
        for k in depths:
            # basic growth estimate: (copies_count^k) * base_size + sum of blockers per level
            for offsets in offsets_set:
                for blockers in blockers_templates:
                    for translation in translations:
                        for anchor in blocker_anchors:
                            for extra in extra_copies_opts:
                                copies = len(offsets) + (extra if extra else 0)
                                # very rough estimate (upper bound)
                                est = (copies ** k) * base_size + 4 * (copies ** (k - 1)) * k
                                if est > MAX_INTERVALS:
                                    continue
                                T = build_pattern(
                                    k=k,
                                    base_seed=base,
                                    offsets=offsets,
                                    blockers=blockers,
                                    translation=translation,
                                    blocker_anchor=anchor,
                                    extra_copies=extra,
                                )
                                s, om, cols, n, Tn = evaluate(T)
                                cand = (s, om, cols, n, Tn, T)
                                if best is None:
                                    best = cand
                                else:
                                    # primary: higher score; ties -> fewer intervals -> more colors
                                    if cand[0] > best[0] + 1e-9:
                                        best = cand
                                    elif abs(cand[0] - best[0]) <= 1e-9:
                                        if cand[3] < best[3]:
                                            best = cand
                                        elif cand[3] == best[3] and cand[2] > best[2]:
                                            best = cand

    # Fallback: original baseline if nothing found (shouldn't happen)
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

    # Try to shrink the raw intervals while preserving witness
    raw_best = best[5]
    shrunk = shrink_intervals(raw_best, max_trials=3000)
    final = normalize_intervals(shrunk)
    # Final sanity: ensure omega>0 and non-empty
    if not final:
        return normalize_intervals(best[5])
    return final

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()