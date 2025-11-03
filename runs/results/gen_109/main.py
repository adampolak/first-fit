# EVOLVE-BLOCK-START

from math import gcd
import random

def construct_intervals(iterations=4):
    """
    Build a sequence of open intervals (integer endpoints) aimed at maximizing
    FirstFit/omega by order optimization (greedy adversarial sequencing).
    Returns: list of (l, r) integer pairs (open intervals).
    """
    # ---------------------
    # 1) Build base geometry (Figure-4 style)
    # ---------------------
    def build_figure4(k=4, base_seed=None, offsets=(2, 6, 10, 14),
                      blockers=None, translation='left', blocker_anchor='left'):
        if base_seed is None:
            base_seed = [(0.0, 1.0)]
        if blockers is None:
            blockers = ((1, 5), (12, 16), (4, 9), (8, 13))
        T = list(base_seed)
        for _ in range(k):
            lo = min(l for l, r in T)
            hi = max(r for l, r in T)
            delta = hi - lo
            center = (lo + hi) / 2.0

            S = []
            for start in offsets:
                if translation == 'left':
                    off = delta * start - lo
                else:
                    off = delta * start - center
                for (l, r) in T:
                    S.append((l + off, r + off))

            for (a, b) in blockers:
                if blocker_anchor == 'left':
                    S.append((delta * a, delta * b))
                else:
                    S.append((delta * a - center, delta * b - center))

            T = S
        return T

    # ---------------------
    # 2) Normalize endpoints to compact even integer grid
    # ---------------------
    def normalize_grid(intervals):
        if not intervals:
            return []
        pts = sorted(set(x for seg in intervals for x in seg))
        coord = {}
        cur = 0
        for p in pts:
            coord[p] = cur
            cur += 2
        return [(coord[l], coord[r]) for (l, r) in intervals]

    # ---------------------
    # 3) Overlap/clique/FirstFit simulators (robust)
    # ---------------------
    def overlaps(a, b):
        (l1, r1), (l2, r2) = a, b
        return max(l1, l2) < min(r1, r2)

    def clique_number(intervals):
        events = []
        for (l, r) in intervals:
            if l < r:
                events.append((l, +1))
                events.append((r, -1))
        # open intervals: handle right (-1) before left (+1) at same coordinate
        events.sort(key=lambda e: (e[0], 0 if e[1] == -1 else 1))
        cur = best = 0
        for _, t in events:
            cur += t
            if cur > best:
                best = cur
        return best

    def firstfit_colors(intervals):
        """
        Exact FirstFit: maintain color classes (lists of intervals).
        This is O(n^2) but n is modest for our instances.
        """
        colors = []
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

    # ---------------------
    # 4) Greedy reordering heuristics
    #    At each step pick the remaining interval that would be assigned to the
    #    largest color index under the current partial FirstFit state.
    #    Use per-color 'last_end' fast test (valid because intervals in a color
    #    are pairwise disjoint and appended in arrival order).
    # ---------------------
    def greedy_order(intervals, tie_break='longest', rng=None):
        # intervals: list of (l, r) integer tuples (normalized)
        remaining = list(intervals)
        order = []
        # color_last_end: current right endpoint for each color class
        color_last_end = []

        # helper to compute assigned color index quickly
        def assigned_index(iv):
            l, r = iv
            for i, le in enumerate(color_last_end):
                # open intervals: we can place in color i if iv.l >= last_end
                if l >= le:
                    return i
            return len(color_last_end)

        # pre-shuffle remaining deterministically for tie-breaking stability
        if rng is not None:
            rng.shuffle(remaining)

        while remaining:
            best_idx = -1
            candidates = []
            # scan remaining and compute assigned index under current last_end's
            for iv in remaining:
                idx = assigned_index(iv)
                if idx > best_idx:
                    best_idx = idx
                    candidates = [iv]
                elif idx == best_idx:
                    candidates.append(iv)

            # tie-breaking among candidates
            if len(candidates) == 1:
                pick = candidates[0]
            else:
                if tie_break == 'longest':
                    pick = max(candidates, key=lambda x: (x[1] - x[0], -x[0]))
                elif tie_break == 'shortest':
                    pick = min(candidates, key=lambda x: (x[1] - x[0], x[0]))
                elif tie_break == 'leftmost':
                    pick = min(candidates, key=lambda x: (x[0], x[1]))
                elif tie_break == 'rightmost':
                    pick = max(candidates, key=lambda x: (x[1], -x[0]))
                elif tie_break.startswith('random'):
                    if rng is None:
                        pick = candidates[0]
                    else:
                        pick = rng.choice(candidates)
                else:
                    pick = candidates[0]

            # assign pick to color classes
            assigned = None
            l, r = pick
            for i, le in enumerate(color_last_end):
                if l >= le:
                    assigned = i
                    color_last_end[i] = r
                    break
            if assigned is None:
                color_last_end.append(r)
            order.append(pick)
            remaining.remove(pick)
        return order

    # ---------------------
    # 5) Conservative deterministic pruning pass
    #    Remove intervals one-by-one (prefer long ones) only if the observed
    #    FirstFit/omega ratio is preserved (no drop).
    # ---------------------
    def prune_preserve_ratio(intervals, target_ratio, max_removals=200):
        cur = list(intervals)
        if not cur:
            return cur
        base_om = clique_number(cur)
        if base_om == 0:
            return cur
        base_cols = firstfit_colors(cur)
        base_ratio = base_cols / base_om
        # acceptance threshold
        threshold = target_ratio if target_ratio is not None else base_ratio
        # attempt removals sorted by descending length
        def length(iv): return iv[1] - iv[0]
        # deterministic ordering stable tie-break
        idx_order = sorted(range(len(cur)), key=lambda i: (-length(cur[i]), cur[i][0], cur[i][1]))
        removals = 0
        changed = True
        while changed and removals < max_removals:
            changed = False
            for i in list(idx_order):
                if i >= len(cur):
                    continue
                cand = cur[:i] + cur[i+1:]
                if not cand:
                    continue
                om = clique_number(cand)
                if om == 0:
                    continue
                cols = firstfit_colors(cand)
                ratio = cols / om
                # accept only if ratio not dropping below threshold
                if ratio >= threshold - 1e-12:
                    cur = cand
                    # recompute idx_order
                    idx_order = sorted(range(len(cur)), key=lambda j: (-length(cur[j]), cur[j][0], cur[j][1]))
                    changed = True
                    removals += 1
                    break
        return cur

    # ---------------------
    # 6) Build, search heuristics, pick best ordering
    # ---------------------
    # Build the canonical figure-4 4-level geometry (this matches prior best
    # geometry that produced 13 colors with omega 5). We keep k=4 by default
    k = 4 if iterations <= 4 else iterations
    raw = build_figure4(k=k, base_seed=[(0.0, 1.0)],
                        offsets=(2, 6, 10, 14),
                        blockers=((1, 5), (12, 16), (4, 9), (8, 13)),
                        translation='left', blocker_anchor='left')

    # normalize to integer grid first (this is the geometry that clique/FF will use)
    base_norm = normalize_grid(raw)
    if not base_norm:
        return []

    base_omega = clique_number(base_norm)
    if base_omega <= 0:
        return base_norm

    # candidate heuristics to try (deterministic seeds)
    heuristics = [
        ('longest', None),
        ('shortest', None),
        ('leftmost', None),
        ('rightmost', None),
        ('random-42', random.Random(42)),
        ('random-101', random.Random(101)),
        ('random-7', random.Random(7)),
    ]

    best_order = None
    best_cols = -1
    best_ratio = -1.0
    best_n = None

    # Always include the original (baseline) order as a candidate
    baseline_cols = firstfit_colors(base_norm)
    baseline_ratio = baseline_cols / base_omega
    best_order = list(base_norm)
    best_cols = baseline_cols
    best_ratio = baseline_ratio
    best_n = len(base_norm)

    # try each greedy heuristic
    for (name, rng) in heuristics:
        tb = name.split('-')[0]  # tie_break key
        order = greedy_order(base_norm, tie_break=tb, rng=rng)
        cols = firstfit_colors(order)
        ratio = cols / base_omega
        n = len(order)
        # choose better by ratio; tie-break fewer intervals then more colors
        if ratio > best_ratio + 1e-12 or (abs(ratio - best_ratio) <= 1e-12 and (n < best_n or (n == best_n and cols > best_cols))):
            best_order = order
            best_cols = cols
            best_ratio = ratio
            best_n = n

    # Conservative pruning: only if it preserves ratio (do not drop)
    pruned = prune_preserve_ratio(best_order, target_ratio=best_ratio, max_removals=200)
    # re-evaluate after pruning
    pruned_cols = firstfit_colors(pruned)
    pruned_om = clique_number(pruned)
    pruned_ratio = pruned_cols / pruned_om if pruned_om > 0 else -1.0

    # Accept pruned if ratio not decreased (minor numeric tolerance)
    if pruned_ratio >= best_ratio - 1e-12:
        final = pruned
    else:
        final = best_order

    # Final normalization (endpoints already integer grid); ensure tuple ints
    final_int = [(int(l), int(r)) for (l, r) in final]
    return final_int

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()