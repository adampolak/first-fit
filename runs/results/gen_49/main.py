# EVOLVE-BLOCK-START

from math import gcd

# ------------------------------
# Geometry and overlap utilities
# ------------------------------

def overlaps(a, b):
    """Open-interval overlap test: return True iff intervals overlap."""
    (l1, r1), (l2, r2) = a, b
    return max(l1, l2) < min(r1, r2)

def clique_number(intervals):
    """
    Compute omega (maximum number of intervals covering a single point) using sweep.
    Open intervals: treat right endpoints before left endpoints.
    """
    events = []
    for (l, r) in intervals:
        if l >= r:
            continue
        events.append((l, +1))
        events.append((r, -1))
    # sort with right(-1) before left(+1) at ties
    events.sort(key=lambda e: (e[0], 0 if e[1] == -1 else 1))
    cur = 0
    best = 0
    prev_x = None
    for x, t in events:
        cur += t
        if cur > best:
            best = cur
        prev_x = x
    return best

def coverage_cells(intervals):
    """
    Return list of coverage 'cells': [(l_i, r_i, m_i)], where for any x in (l_i, r_i),
    exactly m_i intervals cover x. Open endpoints: cells span between consecutive endpoints.
    """
    events = []
    for (l, r) in intervals:
        if l >= r:
            continue
        events.append((l, +1))
        events.append((r, -1))
    if not events:
        return []
    # sort with right(-1) before left(+1) to honor open intervals
    events.sort(key=lambda e: (e[0], 0 if e[1] == -1 else 1))
    # unique sorted endpoints
    xs = sorted(set(x for x, _ in events))
    # sweep to compute coverage level just after each endpoint
    # We compute coverage value on each open cell (xs[i], xs[i+1]).
    idx_events = 0
    cur = 0
    cells = []
    # Preprocess: for each endpoint, accumulate delta before proceeding to next coordinate
    # For open-interval semantics we apply all -1 first then +1 at each coordinate.
    # Our event ordering already ensures that.
    # Compute coverage between successive xs.
    counts_at = {}
    i = 0
    while i < len(xs):
        x = xs[i]
        while idx_events < len(events) and events[idx_events][0] == x:
            cur += events[idx_events][1]
            idx_events += 1
        counts_at[x] = cur
        i += 1
    for i in range(len(xs) - 1):
        a, b = xs[i], xs[i + 1]
        # coverage on any interior point equals counts_at at right after processing at 'a'
        m = counts_at[a]
        if a < b:
            cells.append((a, b, m))
    return cells

def extract_runs(cells, cap):
    """
    From coverage cells, extract maximal runs where coverage <= cap.
    Returns list of (L, R, min_m_on_run, max_m_on_run).
    """
    runs = []
    cur_L = None
    cur_R = None
    cur_min = None
    cur_max = None
    for (l, r, m) in cells:
        if m <= cap:
            if cur_L is None:
                cur_L, cur_R = l, r
                cur_min = m
                cur_max = m
            else:
                # extend contiguous if adjacent
                if abs(l - cur_R) < 1e-12:
                    cur_R = r
                    if m < cur_min:
                        cur_min = m
                    if m > cur_max:
                        cur_max = m
                else:
                    runs.append((cur_L, cur_R, cur_min, cur_max))
                    cur_L, cur_R = l, r
                    cur_min = m
                    cur_max = m
        else:
            if cur_L is not None:
                runs.append((cur_L, cur_R, cur_min, cur_max))
                cur_L = cur_R = cur_min = cur_max = None
    if cur_L is not None:
        runs.append((cur_L, cur_R, cur_min, cur_max))
    return runs

# ------------------------------
# FirstFit partition and helpers
# ------------------------------

def firstfit_partition(intervals):
    """
    Simulate FirstFit, return color classes in arrival order.
    colors: list of lists (intervals).
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
    return colors

def firstfit_colors(intervals):
    return len(firstfit_partition(intervals))

# ------------------------------
# Normalization
# ------------------------------

def normalize_intervals(intervals):
    """
    Normalize to small integer coordinates:
    - scale by 2 and round (to clear 0.5),
    - translate so min coord is 0,
    - divide by gcd of endpoints to shrink if possible.
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

# ------------------------------
# Baseline builder (4-copy + 4-blocker)
# ------------------------------

def make_copies(T, offsets, delta, lo):
    S = []
    for start in offsets:
        offset = delta * start - lo
        for (l, r) in T:
            S.append((l + offset, r + offset))
    return S

def build_baseline(k=4, offsets=(2, 6, 10, 14), blockers=((1, 5), (12, 16), (4, 9), (8, 13))):
    """
    Canonical 4-copy/4-blocker recursion (Figure 4 style).
    """
    T = [(0.0, 1.0)]
    for _ in range(k):
        lo = min(l for l, r in T)
        hi = max(r for l, r in T)
        delta = hi - lo
        S = []
        S.extend(make_copies(T, offsets, delta, lo))
        for (a, b) in blockers:
            S.append((delta * a, delta * b))
        T = S
    return normalize_intervals(T)

# ------------------------------
# Corridor-waves augmenter
# ------------------------------

def colors_hit_by_segment(segment, color_classes):
    """
    Return number of color classes that have at least one interval overlapping the segment.
    """
    (L, R) = segment
    hit = 0
    for cls in color_classes:
        ok = False
        for (l, r) in cls:
            if max(l, L) < min(r, R):
                ok = True
                break
        if ok:
            hit += 1
    return hit

def pick_augmentable_run(baseline, omega_cap):
    """
    Given a baseline instance and omega, find a run R where:
      - coverage <= omega_cap-1 (so adding a single wave keeps omega <= omega_cap),
      - it intersects every color class (so a wave placed here forces a new color).
    Additionally, prefer runs whose interior contains a cell with coverage <= omega_cap-2,
    enabling two overlapping waves (sharing a small core).
    Returns (run_L, run_R, allow_two_waves, core_cell) or None.
    """
    cells = coverage_cells(baseline)
    # runs where baseline coverage <= omega-1
    runs = extract_runs(cells, omega_cap - 1)
    color_classes = firstfit_partition(baseline)
    total_colors = len(color_classes)

    best = None
    # Score runs by (covers_all_colors, allow_two_waves, length, hit_count)
    for (L, R, min_m, max_m) in runs:
        if R - L <= 0:
            continue
        hit = colors_hit_by_segment((L, R), color_classes)
        covers_all = (hit == total_colors)
        # find an internal cell with m <= omega-2 to allow 2-wave overlap core
        allow_two = False
        core = None
        if min_m <= omega_cap - 2:
            # pick the widest cell within (L,R) that has m <= omega-2
            best_len = -1.0
            for (a, b, m) in cells:
                if m <= omega_cap - 2:
                    # ensure (a,b) lies within (L,R)
                    lo = max(L, a)
                    hi = min(R, b)
                    if hi > lo:
                        length = hi - lo
                        if length > best_len:
                            best_len = length
                            # keep a small margin inside the cell to avoid touching endpoints
                            eps = max(1.0, (hi - lo) / 10.0)
                            core = (lo + eps/2.0, hi - eps/2.0) if (hi - lo) > 2*eps else (lo + (hi-lo)/4.0, hi - (hi-lo)/4.0)
            allow_two = core is not None
        score = (
            1 if covers_all else 0,
            1 if allow_two else 0,
            R - L,
            hit
        )
        if best is None or score > best[0]:
            best = (score, (L, R, allow_two, core, hit, total_colors))
    if best is None:
        return None
    _, (L, R, allow_two, core, hit, total) = best
    if hit < total:
        return None
    # ensure we have some core; if not allowing two waves, fabricate a tiny core
    if not allow_two:
        # build a tiny core inside (L,R) where coverage <= omega-1 (always true by construction)
        mid = (L + R) / 2.0
        core = (mid - (R - L) / 100.0, mid + (R - L) / 100.0)
    return (L, R, allow_two, core)

def add_wave(intervals, segment, omega_cap):
    """
    Try to append 'segment' as a wave; accept only if omega remains <= omega_cap
    and it forces a new color (overlaps one interval from every existing color class).
    Returns (new_intervals, success).
    """
    cur_colors = firstfit_partition(intervals)
    total_colors = len(cur_colors)
    # must hit all colors
    if colors_hit_by_segment(segment, cur_colors) < total_colors:
        return intervals, False
    cand = intervals + [segment]
    if clique_number(cand) > omega_cap:
        return intervals, False
    # also check FirstFit assigns a new color
    if firstfit_colors(cand) <= total_colors:
        return intervals, False
    return cand, True

def augment_with_corridor_waves(baseline):
    """
    Main augmenter:
    - Find an augmentable run R intersecting all colors,
    - Add one wave J1 = (L, c_hi),
    - If possible, add a second wave J2 = (c_lo, R),
      both sharing only a tiny 'core' segment C to keep omega stable.
    """
    if not baseline:
        return baseline
    omega_cap = clique_number(baseline)
    if omega_cap <= 0:
        return baseline

    picked = pick_augmentable_run(baseline, omega_cap)
    if picked is None:
        return baseline
    (L, R, allow_two, core) = picked

    # Define a small core safely inside (L,R)
    c_lo, c_hi = core
    # Ensure core lies strictly inside (L,R)
    c_lo = max(L + 1.0, min(c_lo, R - 2.0))
    c_hi = min(R - 1.0, max(c_hi, L + 2.0))
    if not (L < c_lo < c_hi < R):
        return baseline

    # Construct first wave taking left wing + core
    J1 = (L, c_hi)
    T1, ok1 = add_wave(baseline, J1, omega_cap)
    if not ok1:
        # try right wing first if left fails
        J1_alt = (c_lo, R)
        T1, ok1 = add_wave(baseline, J1_alt, omega_cap)
        if not ok1:
            return baseline  # no augmentation possible
        else:
            # right-wing succeeded
            J1 = J1_alt

    # Try a second wave: the complementary wing
    J2 = (c_lo, R) if J1 == (L, c_hi) else (L, c_hi)
    # Try strictly limiting overlap to the core by slightly tightening towards core
    tight = 0.0
    T2, ok2 = add_wave(T1, J2, omega_cap)
    if not ok2 and allow_two:
        # progressively shrink to keep overlap minimal while attempting to cover all colors
        for shrink in [0.0, 0.1, 0.2, 0.3]:
            if J1 == (L, c_hi):
                J2_try = (c_lo + shrink * (c_hi - c_lo), R)
            else:
                J2_try = (L, c_hi - shrink * (c_hi - c_lo))
            T2, ok2 = add_wave(T1, J2_try, omega_cap)
            if ok2:
                break

    return T2 if ok2 else T1

# ------------------------------
# Orchestrator
# ------------------------------

def construct_intervals():
    """
    Build a strong baseline using the 4-copy/4-blocker recursion,
    then augment it with up to two corridor waves that raise FirstFit colors
    while preserving the clique number.
    """
    # Candidate baselines: vary offsets slightly to improve chance
    depths = [3, 4, 5]
    offsets_set = [
        (2, 6, 10, 14),
        (1, 5, 9, 13),
        (3, 7, 11, 15),
        (0, 4, 8, 12),
    ]
    blockers = ((1, 5), (12, 16), (4, 9), (8, 13))

    best = None  # (score, omega, cols, n, intervals)
    for k in depths:
        for offs in offsets_set:
            base = build_baseline(k=k, offsets=offs, blockers=blockers)
            # augment
            aug = augment_with_corridor_waves(base)
            # Evaluate
            om = clique_number(aug)
            cols = firstfit_colors(aug)
            ratio = cols / max(1, om)
            score = ratio - 1e-6 * (len(aug) / 10000.0)
            cand = (score, om, cols, len(aug), aug)
            if best is None:
                best = cand
            else:
                if cand[0] > best[0] + 1e-9:
                    best = cand
                elif abs(cand[0] - best[0]) <= 1e-9:
                    # smaller n preferred, then more colors
                    if cand[3] < best[3] or (cand[3] == best[3] and cand[2] > best[2]):
                        best = cand

    # Fallback to baseline with default k if no augmented candidate improves
    return best[4] if best is not None else build_baseline(k=4)

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()