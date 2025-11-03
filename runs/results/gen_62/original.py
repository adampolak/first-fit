# EVOLVE-BLOCK-START

from math import gcd

# ------------------------------
# Overlap and sweep utilities
# ------------------------------

def overlaps(a, b):
    """Open-interval overlap test: return True iff intervals overlap."""
    (l1, r1), (l2, r2) = a, b
    return max(l1, l2) < min(r1, r2)

def clique_number(intervals):
    """
    Compute omega (maximum number of intervals covering a single point) using sweep.
    Open intervals: endpoints do not contribute; handle right(-1) before left(+1) at ties.
    """
    events = []
    for (l, r) in intervals:
        if l >= r:
            continue
        events.append((l, +1))
        events.append((r, -1))
    if not events:
        return 0
    events.sort(key=lambda e: (e[0], 0 if e[1] == -1 else 1))
    cur = 0
    best = 0
    for _, t in events:
        cur += t
        if cur > best:
            best = cur
    return best

# ------------------------------
# FirstFit simulation
# ------------------------------

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

# ------------------------------
# Pattern builder (4-copy + 4-blockers)
# ------------------------------

def make_copies(T, offsets, delta, lo, center, translation):
    """
    Create translated copies of T according to offsets and translation rule.
    translation in {'left', 'center'}.
    """
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
    """
    Add blockers scaled by delta. Anchor may be 'left' (absolute) or 'center' (center-shifted).
    """
    for (a, b) in blockers:
        if anchor == 'left':
            S.append((delta * a, delta * b))
        else:
            S.append((delta * a - center, delta * b - center))
    return S

def build_pattern(k, base_seed, offsets, blockers, translation, blocker_anchor):
    """
    Recursively expand the base_seed k times using the 4-copy + 4-blocker scheme.
    """
    T = list(base_seed)
    for _ in range(k):
        lo = min(l for l, r in T)
        hi = max(r for l, r in T)
        delta = hi - lo
        center = (lo + hi) / 2.0

        S = make_copies(T, offsets, delta, lo, center, translation)
        S = add_blockers(S, blockers, delta, blocker_anchor, center)
        T = S
    return T

# ------------------------------
# Normalization and evaluation
# ------------------------------

def normalize_intervals(intervals):
    """
    Normalize to small integer coordinates:
    - scale by 2 to eliminate .5 if produced by center shifts
    - translate so min coordinate is >= 0
    - divide by global gcd to shrink
    """
    if not intervals:
        return intervals
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
    score = ratio - 1e-6 * (n / 10000.0)
    return (score, om, cols, n, Tn)

# ------------------------------
# Pruning: strict then relaxed
# ------------------------------

def prune_strict(intervals, target_cols, target_omega):
    """
    Strict pruning: remove any interval whose removal keeps both FirstFit colors and omega unchanged.
    Deterministic order: try longer intervals first to remove redundant long blockers if any.
    """
    cur = list(intervals)
    def length(iv): return iv[1] - iv[0]
    order = sorted(range(len(cur)), key=lambda i: (-length(cur[i]), i))
    changed = True
    while changed:
        changed = False
        for idx in list(order):
            if idx >= len(cur):
                continue
            cand = cur[:idx] + cur[idx+1:]
            Tn = normalize_intervals(cand)
            om = clique_number(Tn)
            cols = firstfit_colors(Tn)
            if om == target_omega and cols == target_cols:
                cur = cand
                order = sorted(range(len(cur)), key=lambda i: (-length(cur[i]), i))
                changed = True
                break
    return cur

def prune_relaxed(intervals, min_ratio, max_omega):
    """
    Relaxed pruning: remove intervals if the ratio remains >= min_ratio and omega <= max_omega.
    Prefer to remove longer intervals first; greedily iterate to fixpoint.
    """
    cur = list(intervals)
    def length(iv): return iv[1] - iv[0]
    order = sorted(range(len(cur)), key=lambda i: (-length(cur[i]), i))
    changed = True
    while changed:
        changed = False
        for idx in list(order):
            if idx >= len(cur):
                continue
            cand = cur[:idx] + cur[idx+1:]
            Tn = normalize_intervals(cand)
            om = clique_number(Tn)
            if om == 0 or om > max_omega:
                continue
            cols = firstfit_colors(Tn)
            ratio = cols / om
            if ratio + 1e-12 >= min_ratio and om <= max_omega:
                cur = cand
                order = sorted(range(len(cur)), key=lambda i: (-length(cur[i]), i))
                changed = True
                break
    return cur

# ------------------------------
# Orchestrator
# ------------------------------

def construct_intervals():
    """
    Build a sequence of open intervals that aims to maximize FirstFit/OPT.
    We sweep several recommended blueprints and pick the best validated candidate.
    Then run a conservative two-phase pruning to reduce interval count without hurting the ratio.
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
        [(0.0, 1.0), (2.0, 3.0)],          # two disjoint seeds variant
    ]
    extra_copies_opts = [0, 1, 2]  # allow up to two extra translated copies on first level

    best = None  # (score, om, cols, n, intervals_norm, raw_intervals)
    # Enumerate combinations with a guard on worst-case explosion
    for base in base_seeds:
        for k in depths:
            if (4 ** k) * (len(base) + 2) > 2000:
                continue
            for offsets in offsets_set:
                for blockers in blockers_templates:
                    for translation in translations:
                        for anchor in blocker_anchors:
                            for extra in extra_copies_opts:
                                offs = list(offsets)
                                if extra >= 1:
                                    last = offs[-1] if offs else 0
                                    offs.append(last + 4)
                                if extra == 2:
                                    last = offs[-2] if len(offs) >= 2 else (offs[-1] if offs else 0)
                                    offs.append(last + 8)
                                offs = tuple(offs)

                                T = build_pattern(
                                    k=k,
                                    base_seed=base,
                                    offsets=offs,
                                    blockers=blockers,
                                    translation=translation,
                                    blocker_anchor=anchor,
                                )
                                score, om, cols, n, Tn = evaluate(T)
                                cand = (score, om, cols, n, Tn, T)
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

    # Fallback to baseline if search didn't produce anything
    if best is None:
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

    # Apply two-phase pruning to the raw best intervals
    raw_best = best[5]
    _, best_om, best_cols, _, _, _ = best
    # Phase 1: strict pruning (keep alg and omega identical)
    phase1 = prune_strict(raw_best, best_cols, best_om)
    # Phase 2: relaxed pruning (keep ratio and omega cap)
    target_ratio = best_cols / best_om if best_om > 0 else 0.0
    phase2 = prune_relaxed(phase1, target_ratio, best_om)

    final = normalize_intervals(phase2)
    if not final:
        return best[4]
    # Final safety: if ratio dropped (shouldn't), revert to best normalized
    om = clique_number(final)
    cols = firstfit_colors(final)
    if om == 0 or cols / om + 1e-12 < target_ratio or om > best_om:
        return best[4]
    return final

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()