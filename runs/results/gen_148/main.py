# EVOLVE-BLOCK-START

from math import gcd

def construct_intervals(iterations=4):
    """
    Build an adversarial interval sequence for FirstFit.
    Parameter sweep over depths, offsets rotation, blocker templates,
    translations, blocker anchors, and normalization cadence.
    Post-process with conservative shrink_prune.
    Returns a list of (l,r) with integer endpoints.
    """

    # --- Core utilities ---

    def firstfit_colors(intervals):
        """Simulate FirstFit by tracking last endpoint per color."""
        last_end = []
        for l, r in intervals:
            placed = False
            for i, le in enumerate(last_end):
                if l >= le:
                    last_end[i] = r
                    placed = True
                    break
            if not placed:
                last_end.append(r)
        return len(last_end)

    def clique_number(intervals):
        """Sweep-line to compute ω for open intervals."""
        events = []
        for l, r in intervals:
            if l < r:
                events.append((l, +1))
                events.append((r, -1))
        # process -1 before +1 at ties
        events.sort(key=lambda e: (e[0], 0 if e[1] == -1 else 1))
        cur = best = 0
        for _, t in events:
            cur += t
            if cur > best:
                best = cur
        return best

    def normalize_intervals(intervals):
        """
        Normalize floats to integers for ratio tests:
        scale by 2, round, shift min to 0, divide by gcd.
        """
        if not intervals:
            return []
        scaled = []
        for l, r in intervals:
            L = int(round(l * 2))
            R = int(round(r * 2))
            scaled.append((L, R))
        # shift
        minc = min(min(l, r) for l, r in scaled)
        shifted = [(l - minc, r - minc) for l, r in scaled]
        # divide gcd
        g = 0
        for l, r in shifted:
            g = gcd(g, abs(l))
            g = gcd(g, abs(r))
        if g > 1:
            return [(l // g, r // g) for l, r in shifted]
        else:
            return shifted

    def normalize_grid(intervals):
        """
        Final mapping: unique endpoints -> even integers 0,2,4,...
        """
        if not intervals:
            return []
        pts = sorted({x for seg in intervals for x in seg})
        coord = {}
        cur = 0
        for x in pts:
            coord[x] = cur
            cur += 2
        return [(coord[l], coord[r]) for l, r in intervals]

    def shrink_prune(raw, target_ratio):
        """
        Remove long intervals one by one if FirstFit/ω ≥ target_ratio.
        Deterministic: tries longer intervals first.
        """
        cur = list(raw)
        if not cur:
            return cur
        base_norm = normalize_intervals(cur)
        base_cols = firstfit_colors(base_norm)
        base_om = clique_number(base_norm)
        if base_om == 0:
            return cur
        # target_ratio default to observed
        if target_ratio is None:
            target_ratio = base_cols / base_om

        def length(iv):
            return iv[1] - iv[0]

        # order by decreasing length
        order = sorted(range(len(cur)), key=lambda i: (-length(cur[i]), i))
        changed = True
        while changed:
            changed = False
            for idx in order:
                if idx >= len(cur):
                    continue
                cand = cur[:idx] + cur[idx+1:]
                norm = normalize_intervals(cand)
                if not norm:
                    continue
                cols = firstfit_colors(norm)
                om = clique_number(norm)
                if om > 0 and cols/om >= target_ratio:
                    cur = cand
                    # recompute order
                    order = sorted(range(len(cur)), key=lambda i: (-length(cur[i]), i))
                    changed = True
                    break
        return cur

    # --- Pattern builders ---

    def make_copies(T, offsets, delta, lo, center, translation):
        S = []
        for s in offsets:
            if translation == 'left':
                off = delta * s - lo
            else:
                off = delta * s - center
            for l, r in T:
                S.append((l + off, r + off))
        return S

    def add_blockers(T, blockers, delta, anchor, center):
        S = list(T)
        for a, b in blockers:
            if anchor == 'left':
                S.append((delta * a, delta * b))
            else:
                S.append((delta * a - center, delta * b - center))
        return S

    def build_pattern(k, offsets, blockers, translation, anchor, norm_schedule):
        T = [(0.0, 1.0)]
        for i in range(k):
            lo = min(l for l, r in T)
            hi = max(r for l, r in T)
            delta = hi - lo
            center = (lo + hi) / 2.0
            # optionally normalize per level
            if norm_schedule == 'per_level':
                # scale endpoints to [0,1]
                T = [((l - lo)/delta, (r - lo)/delta) for l, r in T]
                lo, hi = 0.0, 1.0
                delta = 1.0
                center = 0.5
            S = make_copies(T, offsets, delta, lo, center, translation)
            S = add_blockers(S, blockers, delta, anchor, center)
            T = S
        if norm_schedule == 'final':
            # reset span to [0,1]
            lo = min(l for l, r in T)
            hi = max(r for l, r in T)
            span = hi - lo
            if span > 0:
                T = [((l - lo)/span, (r - lo)/span) for l, r in T]
        return T

    # --- Parameter sets ---

    offsets_cycle = [
        (2,6,10,14),  # A
        (1,5,9,13),   # B
        (3,7,11,15),  # C
        (0,4,8,12),   # D
    ]
    blocker_templates = [
        ((1,5),(12,16),(4,9),(8,13)),  # A
        ((0,4),(6,10),(8,12),(14,18)), # B
        ((2,6),(4,8),(10,14),(12,16)), # C
    ]
    translations = ['left','center']
    anchors = ['left','center']
    depths = [3,4,5]
    norm_schedules = ['none','per_level','final']

    best_raw = None
    best_ratio = -1.0
    best_size = 0
    best_cols = best_om = 0

    # --- Search loop ---

    for k in depths:
        # rotate offsets by depth
        offsets = offsets_cycle[(k-3) % len(offsets_cycle)]
        for blockers in blocker_templates:
            for translation in translations:
                for anchor in anchors:
                    for norm_sc in norm_schedules:
                        T = build_pattern(k, offsets, blockers, translation, anchor, norm_sc)
                        n = len(T)
                        if n > 2000:
                            continue
                        om = clique_number(T)
                        if om == 0:
                            continue
                        cols = firstfit_colors(T)
                        ratio = cols / om
                        # select by ratio, then smaller n, then larger cols
                        better = False
                        if ratio > best_ratio + 1e-12:
                            better = True
                        elif abs(ratio - best_ratio) <= 1e-12:
                            if n < best_size or (n == best_size and cols > best_cols):
                                better = True
                        if better:
                            best_ratio, best_raw, best_size, best_cols, best_om = \
                                ratio, T, n, cols, om

    # fallback to classic if none found
    if best_raw is None:
        # standard 4‐iteration Figure‐4
        best_raw = [(0.0,1.0)]
        for _ in range(iterations):
            lo = min(l for l, r in best_raw)
            hi = max(r for l, r in best_raw)
            delta = hi - lo
            S = []
            for s in (2,6,10,14):
                S += [(l+delta*s-lo, r+delta*s-lo) for l,r in best_raw]
            S += [(delta*1,delta*5),(delta*12,delta*16),
                  (delta*4,delta*9),(delta*8,delta*13)]
            best_raw = S

    # prune and normalize to integer grid
    target = best_cols / best_om if best_om > 0 else None
    pruned = shrink_prune(best_raw, target)
    final = normalize_grid(pruned)

    return final

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()