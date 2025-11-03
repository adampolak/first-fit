# EVOLVE-BLOCK-START

from math import gcd

def _normalize_grid(intervals):
    """
    Normalize endpoints to a compact integer grid while preserving order.
    Map each unique endpoint to an increasing even integer (0,2,4,...),
    translate so minimum is 0, and divide by gcd to shrink further.
    """
    if not intervals:
        return []
    pts = sorted({x for seg in intervals for x in seg})
    coord = {}
    cur = 0
    for p in pts:
        coord[p] = cur
        cur += 2
    mapped = [(coord[l], coord[r]) for (l, r) in intervals]
    # shift min to 0
    mn = min(min(a, b) for a, b in mapped)
    shifted = [(a - mn, b - mn) for a, b in mapped]
    # divide by gcd of all coordinates
    g = 0
    for a, b in shifted:
        g = gcd(g, a)
        g = gcd(g, b)
    if g > 1:
        shrunk = [(a // g, b // g) for a, b in shifted]
        return shrunk
    return shifted

def ff_count(intervals):
    """
    Efficient FirstFit coloring for open intervals in given arrival order.
    Maintain last_end per color and place an interval in the first color
    whose last_end <= l (open intervals: equality allowed).
    Returns number of colors used.
    """
    last_end = []  # last end placed in each color
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

def omega_open(intervals):
    """
    Compute clique number (maximum number of open intervals covering a point)
    via sweep-line. For open intervals, right endpoints processed before left
    endpoints at the same coordinate (handled because -1 < +1).
    """
    events = []
    for (l, r) in intervals:
        if l < r:
            events.append((l, +1))
            events.append((r, -1))
    if not events:
        return 0
    events.sort(key=lambda e: (e[0], e[1]))  # -1 before +1 at same x
    cur = best = 0
    for _, t in events:
        cur += t
        if cur > best:
            best = cur
    return best

def _make_copies(T, offsets, delta, lo, center, translation):
    S = []
    for start in offsets:
        if translation == 'left':
            off = delta * start - lo
        else:
            off = delta * start - center
        for (l, r) in T:
            S.append((l + off, r + off))
    return S

def _add_blockers(S, blockers, delta, anchor, center):
    for (a, b) in blockers:
        if anchor == 'left':
            S.append((delta * a, delta * b))
        else:
            S.append((delta * a - center, delta * b - center))
    return S

def build_candidate(depth, base_seed, offsets, blockers, translation='left', blocker_anchor='left', extra_first=False, extra_last=False):
    """
    Build pattern by recursively expanding base_seed 'depth' times.
    Options:
      - offsets: tuple/list of offsets (multipliers) for copies at each level
      - blockers: list/tuple of 4 pairs (a,b)
      - translation: 'left' or 'center' anchoring copies
      - blocker_anchor: 'left' or 'center' anchoring blockers
      - extra_first/extra_last: if True, add an extra copy (start=18) on first/last level
    """
    T = list(base_seed)
    for level in range(depth):
        lo = min(l for l, r in T)
        hi = max(r for l, r in T)
        delta = hi - lo
        center = (lo + hi) / 2.0
        offs = tuple(offsets)
        if extra_first and level == 0:
            if 18 not in offs:
                offs = offs + (18,)
        if extra_last and level == depth - 1:
            if 18 not in offs:
                offs = offs + (18,)
        S = _make_copies(T, offs, delta, lo, center, translation)
        S = _add_blockers(S, blockers, delta, blocker_anchor, center)
        T = S
    return T

def shrink_prune(intervals, target_ratio=None, max_attempts=2000):
    """
    Conservative pruning: try to remove intervals (prefer long ones) while keeping
    FirstFit/omega >= target_ratio (or not decreasing observed ratio if None).
    Works on float endpoints then normalizes at the end.
    Deterministic and stops when no removal accepted or attempts exhausted.
    """
    cur = list(intervals)
    if not cur:
        return cur
    # Evaluate baseline on normalized form
    cur_norm = _normalize_grid(cur)
    base_om = omega_open(cur_norm)
    if base_om == 0:
        return cur
    base_cols = ff_count(cur_norm)
    base_ratio = base_cols / base_om
    if target_ratio is None:
        target_ratio = base_ratio
    # Order candidate removals by descending length (prefer removing long blockers)
    def length(iv): return iv[1] - iv[0]
    order = sorted(range(len(cur)), key=lambda i: (-length(cur[i]), i))
    attempts = 0
    changed = True
    while changed and attempts < max_attempts:
        changed = False
        attempts += 1
        for idx in list(order):
            if idx >= len(cur):
                continue
            cand = cur[:idx] + cur[idx+1:]
            cand_norm = _normalize_grid(cand)
            om = omega_open(cand_norm)
            if om == 0:
                continue
            cols = ff_count(cand_norm)
            ratio = cols / om
            # Accept removal only if ratio >= target_ratio (no degradation)
            if ratio >= target_ratio - 1e-12:
                cur = cand
                # recompute order and restart scanning
                order = sorted(range(len(cur)), key=lambda i: (-length(cur[i]), i))
                changed = True
                break
    return cur

def construct_intervals():
    """
    Master constructor: search a small parameter space of Figure-4 style expansions,
    evaluate candidates (normalize->omega/FF) and pick the best one. Then apply
    conservative pruning on the raw float representation to reduce redundant blockers.
    """
    # Base seeds
    base_seeds = [
        [(0.0, 1.0)],
        [(0.0, 1.0), (2.0, 3.0)],
    ]
    offsets_candidates = [
        (2, 6, 10, 14),
        (1, 5, 9, 13),
        (3, 7, 11, 15),
        (0, 4, 8, 12),
    ]
    blockers_templates = [
        ((1, 5), (12, 16), (4, 9), (8, 13)),      # canonical
        ((1, 5), (11, 15), (4, 9), (8, 13)),      # slight variation
        ((2, 6), (12, 16), (4, 9), (8, 13)),      # shift first
    ]
    translations = ['left', 'center']
    blocker_anchors = ['left', 'center']
    depths = [3, 4]  # keep search compact and practical

    best = None  # (ratio_score, om, cols, n, normalized_intervals, raw_intervals)
    # Enumerate but guard explosion
    for base in base_seeds:
        for depth in depths:
            # very rough size cap: don't expand if expected size huge
            if (4 ** depth) * (len(base) + 2) > 5000:
                continue
            for offsets in offsets_candidates:
                for blockers in blockers_templates:
                    for translation in translations:
                        for anchor in blocker_anchors:
                            # try with/without an extra first-level copy
                            for extra_first in (False, True):
                                raw = build_candidate(depth=depth, base_seed=base, offsets=offsets, blockers=blockers, translation=translation, blocker_anchor=anchor, extra_first=extra_first)
                                if not raw:
                                    continue
                                norm = _normalize_grid(raw)
                                om = omega_open(norm)
                                if om == 0:
                                    continue
                                cols = ff_count(norm)
                                ratio = cols / om
                                score = ratio  # primary objective
                                cand = (score, om, cols, len(norm), norm, raw)
                                if best is None:
                                    best = cand
                                else:
                                    # prefer higher score, then fewer intervals, then more colors
                                    if cand[0] > best[0] + 1e-12:
                                        best = cand
                                    elif abs(cand[0] - best[0]) <= 1e-12:
                                        if cand[3] < best[3]:
                                            best = cand
                                        elif cand[3] == best[3] and cand[2] > best[2]:
                                            best = cand

    # Fallback to canonical single-depth construction if nothing found
    if best is None:
        T = [(0.0, 1.0)]
        for _ in range(4):
            lo = min(l for l, r in T)
            hi = max(r for l, r in T)
            delta = hi - lo
            S = []
            for start in (2, 6, 10, 14):
                off = delta * start - lo
                for (l, r) in T:
                    S.append((l + off, r + off))
            S += [
                (delta * 1,  delta * 5),
                (delta * 12, delta * 16),
                (delta * 4,  delta * 9),
                (delta * 8,  delta * 13)
            ]
            T = S
        return _normalize_grid(T)

    # Apply conservative pruning on the raw float best to remove redundancies
    best_score, best_om, best_cols, best_n, best_norm, best_raw = best
    target_ratio = best_cols / best_om if best_om > 0 else None
    pruned_raw = shrink_prune(best_raw, target_ratio=target_ratio)
    final = _normalize_grid(pruned_raw)
    # If pruning removed everything or made things worse, fall back to measured best normalized
    if not final:
        return best_norm
    # Final sanity check: ensure we didn't reduce the ratio
    final_om = omega_open(final)
    if final_om == 0:
        return best_norm
    final_cols = ff_count(final)
    if final_cols / final_om + 1e-12 < best_score:
        return best_norm
    return final

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()