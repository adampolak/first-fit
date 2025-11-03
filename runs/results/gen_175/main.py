# EVOLVE-BLOCK-START

from math import gcd
from functools import reduce

# --- Utilities: normalization, FirstFit, omega (open intervals) ---

def _normalize_grid(intervals):
    """
    Normalize endpoints to a compact integer grid while preserving order.
    Each unique endpoint is mapped to an increasing even integer.
    Returns a new list of (l, r) with integer coordinates.
    """
    endpoints = sorted(set([x for seg in intervals for x in seg]))
    coord = {}
    cur = 0
    for e in endpoints:
        coord[e] = cur
        cur += 2  # spacing by 2 to keep even gaps
    return [(coord[l], coord[r]) for (l, r) in intervals]


def _firstfit_count(intervals):
    """
    Fast FirstFit color count using per-color last end tracking.
    For open intervals, a color c is compatible if l >= last_end[c].
    intervals are processed in arrival order.
    """
    last_end = []
    for (l, r) in intervals:
        placed = False
        for i, le in enumerate(last_end):
            if l >= le:
                last_end[i] = r
                placed = True
                break
        if not placed:
            last_end.append(r)
    return len(last_end)


def _omega_open(intervals):
    """
    Compute omega (max number of intervals covering a single point) for open intervals
    via sweep. Process right endpoints before left endpoints on ties.
    """
    events = []
    for (l, r) in intervals:
        if l < r:
            events.append((l, +1))
            events.append((r, -1))
    # right (-1) before left (+1) at the same coordinate
    events.sort(key=lambda e: (e[0], 0 if e[1] == -1 else 1))
    cur = best = 0
    for _, t in events:
        cur += t
        if cur > best:
            best = cur
    return best


# --- Generator primitives: copies, blockers, waves ---

# canonical blocker template (multipliers of delta)
_CANON_BLOCKERS = [(1, 5), (12, 16), (4, 9), (8, 13)]

# Four deterministic start patterns (cycle among these by level)
_PATTERN_A = (2, 6, 10, 14)
_PATTERN_B = (1, 5, 9, 13)
_PATTERN_C = (3, 7, 11, 15)
_PATTERN_D = (0, 4, 8, 12)
_PATTERN_CYCLE = [_PATTERN_A, _PATTERN_B, _PATTERN_C, _PATTERN_D]

def _make_copies(T, starts, delta, lo, center, anchor='left', eps=0.0):
    """
    Make translated copies of T according to starts multipliers.
    If anchor == 'left' offsets use -lo; if 'center' use -center.
    eps is a tiny fraction (0..0.3) applied to each offset deterministically to perturb symmetry.
    """
    S = []
    for start in starts:
        if anchor == 'left':
            offset = delta * start - lo + delta * eps
        else:
            offset = delta * start - center + delta * eps
        for (l, r) in T:
            S.append((l + offset, r + offset))
    return S

def _make_blockers(delta, blockers=_CANON_BLOCKERS, center=0.0, anchor='left', eps=0.0):
    """
    Return scaled blockers (list of (l,r)) using anchor and eps perturbation.
    eps shifts blockers slightly (deterministic).
    """
    S = []
    for (a, b) in blockers:
        if anchor == 'left':
            S.append((delta * a + delta * eps, delta * b + delta * eps))
        else:
            S.append((delta * a - center + delta * eps, delta * b - center + delta * eps))
    return S

def _insert_short_wave(S, lo, delta, shift=0.2, length_frac=0.08, count=3):
    """
    Insert a deterministic short wave: few very short intervals placed in one gap
    between the first two copies (so they overlap many early colors but don't create big cliques).
    count: number of short intervals in the wave; they are pairwise disjoint.
    length_frac: fraction of delta used as length for wave intervals (small).
    shift: fraction into the gap where wave starts.
    """
    if count <= 0:
        return []
    gap_start = lo + delta * 2  # between copy at start=2 and its right neighbor in the canonical layout
    segment_len = delta * 2  # use a span inside which to place small waves
    wlen = max(1e-6, delta * length_frac)
    base = gap_start + segment_len * shift
    wave = []
    spacing = wlen * 1.5
    for i in range(count):
        l = base + i * spacing
        r = l + wlen
        wave.append((l, r))
    return wave

def _insert_long_wave(S, lo, delta, start_frac=0.4, length_frac=0.6, count=2):
    """
    Insert longer waves: slightly longer intervals that bridge copies, designed to
    touch many active colors but still avoid piling up at a single point.
    count: typically 1 or 2 long intervals.
    """
    if count <= 0:
        return []
    span_start = lo + delta * 4  # choose a different gap (deterministic)
    wlen = max(1e-6, delta * length_frac)
    base = span_start + delta * start_frac
    wave = []
    spacing = wlen * 1.1
    for i in range(count):
        l = base + i * spacing
        r = l + wlen
        wave.append((l, r))
    return wave

# --- Build recursive template with cycle + waves + tiny perturbations ---

def build_recursive_cycle(base_T=None, depth=4, extra_first=False, translation='left', wave_mode='none', eps_cycle=(0.0, 0.15, 0.25)):
    """
    Build recursive template:
      - cycle among pattern offsets by level using _PATTERN_CYCLE
      - apply tiny deterministic perturbation eps from eps_cycle using level index
      - optionally add extra copy on first level (start=18)
      - wave_mode: 'none' | 'short' | 'long' | 'both'
    """
    if base_T is None:
        base_T = [(0.0, 1.0)]
    T = list(base_T)
    for level in range(depth):
        lo = min(l for l, r in T)
        hi = max(r for l, r in T)
        delta = hi - lo
        center = (lo + hi) / 2.0

        # choose pattern by level cycle
        starts = list(_PATTERN_CYCLE[level % len(_PATTERN_CYCLE)])
        if extra_first and level == 0:
            # append a deterministic extra copy (start 18) to strengthen initial coupling
            if 18 not in starts:
                starts = starts + [18]

        # choose eps deterministically from eps_cycle
        eps = eps_cycle[level % len(eps_cycle)]

        S = _make_copies(T, starts, delta, lo, center, anchor=translation, eps=eps)
        # append blockers scaled/anchored with epsilon perturbation (use opposite eps to diversify)
        blk_eps = eps_cycle[(level + 1) % len(eps_cycle)]
        blockers = _make_blockers(delta, blockers=_CANON_BLOCKERS, center=center, anchor=translation, eps=blk_eps)
        S.extend(blockers)

        # insert deterministic waves if requested (place them after copies but before blockers)
        if wave_mode in ('short', 'both'):
            S.extend(_insert_short_wave(S, lo, delta, shift=0.18, length_frac=0.06, count=3))
        if wave_mode in ('long', 'both'):
            S.extend(_insert_long_wave(S, lo, delta, start_frac=0.35, length_frac=0.5, count=2))

        T = S
    return T


# --- Deterministic two-stage pruning (witness preserving) ---

def _firstfit_assignment(intervals):
    """
    Return (num_colors, assignment_list) with 0-based color for each interval in order.
    Uses the standard FirstFit policy with per-color last_end.
    """
    last_end = []
    assignment = []
    for (l, r) in intervals:
        placed = False
        for c, le in enumerate(last_end):
            if l >= le:
                last_end[c] = r
                assignment.append(c)
                placed = True
                break
        if not placed:
            last_end.append(r)
            assignment.append(len(last_end) - 1)
    return len(last_end), assignment

def shrink_prune_deterministic(intervals, target_ratio):
    """
    Deterministic two-stage pruning:
      Stage 1: greedy remove any interval (preferring long ones) whose removal keeps ratio >= target_ratio.
      Stage 2: compute frontier intervals (those active at times we saw the maximum color assignment)
               and then attempt further removals but never remove frontier intervals.
    """
    cur = list(intervals)

    def length(iv): return iv[1] - iv[0]
    # Stage 1
    changed = True
    while changed:
        changed = False
        order = sorted(range(len(cur)), key=lambda i: (-length(cur[i]), i))
        for idx in order:
            cand = cur[:idx] + cur[idx+1:]
            if not cand:
                continue
            norm = _normalize_grid(cand)
            om = _omega_open(norm)
            if om == 0:
                continue
            cols = _firstfit_count(norm)
            if cols / om >= target_ratio - 1e-12:
                cur = cand
                changed = True
                break

    # Stage 2: identify frontier intervals (those that contribute to peak color positions)
    # We'll compute assignment on normalized cur and track intervals that are "in play" at times when
    # the online algorithm used its maximum number of colors (i.e., those colored with max_color-1 or
    # those overlapping them).
    norm_cur = _normalize_grid(cur)
    total_cols, assignment = _firstfit_assignment(norm_cur)
    if total_cols <= 1:
        return cur  # nothing to prune further meaningfully

    max_color = total_cols - 1
    # indices in normalized array where assigned color == max_color
    frontier_indices = [i for i, c in enumerate(assignment) if c == max_color]
    # map normalized intervals back to original ones via value->index mapping is cumbersome; instead,
    # recompute frontier by computing overlaps in original coordinates with the time windows
    # where the normalized frontier intervals lie.
    # Get the normalized coordinate -> original coordinate mapping via sorted endpoints
    endpoints = sorted(set([x for seg in cur for x in seg]))
    coord_map = {e: i*2 for i, e in enumerate(endpoints)}  # deterministic mapping used by _normalize_grid
    # Build normalized intervals for cur again (we already have norm_cur)
    # For every frontier normalized interval pick its mid-point in normalized coords and find which original intervals cover that point.
    frontier_pts = []
    for idx in frontier_indices:
        nl, nr = norm_cur[idx]
        if nl < nr:
            mid = (nl + nr) / 2.0
            frontier_pts.append(mid)
    # construct a set of original intervals that cover any frontier mid-point
    frontier_set = set()
    # compute inverse mapping from normalized coords back to original coordinates:
    # Build list of unique endpoints sorted; normalized mapping was even integers starting at 0 with step 2
    uniq = sorted(set([x for seg in cur for x in seg]))
    norm_to_orig = {i*2: uniq[i] for i in range(len(uniq))}
    # For each frontier point (in normalized), map to original approximate coordinate
    for p in frontier_pts:
        # p is even integer or half; map to nearest lower even integer then to original
        key = int(round(p))
        # ensure even
        key = key - (key % 2)
        if key in norm_to_orig:
            orig_pt = norm_to_orig[key]
        else:
            # fallback: take average of endpoints
            orig_pt = (uniq[0] + uniq[-1]) / 2.0
        # collect intervals in cur that cover orig_pt (open intervals)
        for i, (l, r) in enumerate(cur):
            if l < orig_pt < r:
                frontier_set.add(i)

    # Now try to remove non-frontier intervals further
    changed = True
    while changed:
        changed = False
        order = sorted(range(len(cur)), key=lambda i: (-length(cur[i]), i))
        for idx in order:
            if idx in frontier_set:
                continue
            cand = cur[:idx] + cur[idx+1:]
            if not cand:
                continue
            norm = _normalize_grid(cand)
            om = _omega_open(norm)
            if om == 0:
                continue
            cols = _firstfit_count(norm)
            if cols / om >= target_ratio - 1e-12:
                # Accepted removal; need to rebuild frontier_set conservatively (do not recompute frontier to avoid churn)
                # Instead, adjust indices in frontier_set >= idx by -1
                new_frontier = set()
                for fi in frontier_set:
                    if fi < idx:
                        new_frontier.add(fi)
                    elif fi > idx:
                        new_frontier.add(fi - 1)
                frontier_set = new_frontier
                cur = cand
                changed = True
                break

    return cur


# --- Memoization cache for (signature → (cols, omega)) to speed search --- 

_evaluate_cache = {}

def _evaluate_signature(raw_intervals):
    """
    Normalize and compute (cols, omega). Cache by signature (tuple of rounded endpoints).
    """
    # signature: tuple of floats quantized to rational-ish representation to keep keys small
    sig = tuple(round(x, 9) for seg in raw_intervals for x in seg)
    if sig in _evaluate_cache:
        return _evaluate_cache[sig]
    norm = _normalize_grid(raw_intervals)
    om = _omega_open(norm)
    cols = _firstfit_count(norm) if om > 0 else 0
    _evaluate_cache[sig] = (cols, om, norm)
    return _evaluate_cache[sig]


# --- High-level search/sweep (deterministic, small) ---

def construct_intervals():
    """
    Deterministic generator that:
      - sweeps small parameter set (depths, extra_first, wave_mode, translation, eps cycle)
      - builds candidate recursive templates using cycle patterns
      - evaluates FirstFit/omega, chooses best (ratio, then smaller n, then larger cols)
      - applies deterministic two-stage pruning to shrink while preserving ratio
      - returns normalized final intervals
    """
    best = None  # tuple: (ratio, cols, om, n, normalized_intervals, raw_intervals, meta)

    depths = [3, 4, 5]  # small sweep
    extra_first_opts = [False, True]
    wave_modes = ['none', 'short', 'long', 'both']
    translations = ['left', 'center']
    eps_cycles = [
        (0.0, 0.15, 0.25),
        (0.0, 0.1, 0.2),
        (0.0, 0.2, 0.3)
    ]

    # Guard to avoid explosion: ensure estimated size <= 3000
    for depth in depths:
        for extra_first in extra_first_opts:
            for wave_mode in wave_modes:
                for translation in translations:
                    for eps_cycle in eps_cycles:
                        # rough size estimate
                        est = (len(_PATTERN_A) ** depth) * (1 + (1 if extra_first else 0))
                        # be conservative
                        if 4 ** depth * 8 > 4000:
                            continue
                        raw = build_recursive_cycle(depth=depth, extra_first=extra_first,
                                                   translation=translation, wave_mode=wave_mode, eps_cycle=eps_cycle)
                        cols, om, norm = _evaluate_signature(raw)
                        if om == 0:
                            continue
                        ratio = cols / om
                        n = len(norm)
                        meta = dict(depth=depth, extra_first=extra_first, wave_mode=wave_mode, translation=translation, eps_cycle=eps_cycle)
                        cand = (ratio, cols, om, n, norm, raw, meta)
                        if best is None:
                            best = cand
                        else:
                            br = best[0]
                            if ratio > br + 1e-12:
                                best = cand
                            elif abs(ratio - br) <= 1e-12:
                                # tie-break: smaller n, then larger cols
                                if n < best[3] or (n == best[3] and cols > best[1]):
                                    best = cand

    if best is None:
        # Fallback to a safe canonical construction (depth=4)
        raw = build_recursive_cycle(depth=4, extra_first=True, translation='left', wave_mode='none')
        cols, om, norm = _evaluate_signature(raw)
        best = (cols / max(1, om), cols, om, len(norm), norm, raw, dict(fallback=True))

    # Apply deterministic pruning preserving measured ratio
    ratio, cols, om, n, norm_intervals, raw_intervals, meta = best
    target_ratio = cols / om if om > 0 else None
    pruned_raw = shrink_prune_deterministic(raw_intervals, target_ratio)

    final_norm = _normalize_grid(pruned_raw)
    # sanity check: if pruning fails (empty), return the original normalized
    if not final_norm:
        return norm_intervals
    return final_norm


# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()