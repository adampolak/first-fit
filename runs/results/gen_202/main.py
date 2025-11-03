# EVOLVE-BLOCK-START

# Deterministic, parametric construction aimed at maximizing FirstFit/omega
# using a 4-copy recursive backbone with connector "blockers", cycle-of-offsets
# per level, deterministic wave augmentation, and frontier-preserving pruning.

from bisect import bisect_left

# -------------------- Utilities and memoization --------------------

def _firstfit_colors(intervals):
    """
    FirstFit on open intervals in the given arrival order.
    Optimized by tracking only the last end of each color class.
    For interval graphs, this is sufficient: an arriving interval (l,r)
    can be placed in the first color whose last_end <= l.
    """
    # Memoize on exact sequence
    if not hasattr(_firstfit_colors, "_cache"):
        _firstfit_colors._cache = {}
    key = tuple(intervals)
    if key in _firstfit_colors._cache:
        return _firstfit_colors._cache[key]

    last_end = []  # sorted list of last ends per color (monotone non-decreasing per color)
    for (l, r) in intervals:
        # find first color i with last_end[i] <= l
        placed = False
        for i in range(len(last_end)):
            if last_end[i] <= l:
                last_end[i] = r
                placed = True
                break
        if not placed:
            last_end.append(r)
    res = len(last_end)
    _firstfit_colors._cache[key] = res
    return res

def _clique_number(intervals):
    """
    Clique number (maximum number of open intervals covering a point),
    using sweep-line; process right endpoints before left on ties.
    """
    if not hasattr(_clique_number, "_cache"):
        _clique_number._cache = {}
    key = tuple(intervals)
    if key in _clique_number._cache:
        return _clique_number._cache[key]

    events = []
    for (l, r) in intervals:
        if l < r:
            events.append((l, +1))
            events.append((r, -1))
    events.sort(key=lambda e: (e[0], 0 if e[1] == -1 else 1))
    cur = best = 0
    for _, d in events:
        cur += d
        if cur > best:
            best = cur
    _clique_number._cache[key] = best
    return best

def _normalize_grid(intervals):
    """
    Compact integer grid: map unique endpoints to increasing even integers.
    """
    if not intervals:
        return []
    xs = sorted({x for seg in intervals for x in seg})
    coord = {}
    v = 0
    for e in xs:
        coord[e] = v
        v += 2
    return [(coord[l], coord[r]) for (l, r) in intervals]

# -------------------- Backbone builder --------------------

# Offset patterns A,B,C,D used cyclically across levels
_OFFSETS_CYCLE = [
    (2, 6, 10, 14),  # A
    (1, 5, 9, 13),   # B
    (3, 7, 11, 15),  # C
    (0, 4, 8, 12),   # D
]

# Blocker templates: canonical and shifted
_BLOCKERS_TEMPLATES = [
    ((1, 5), (12, 16), (4, 9), (8, 13)),     # canonical
    ((0, 4), (11, 15), (3, 8), (7, 12)),     # shifted
]

def _build_backbone(depth, cycle_start=0, schedule='split', extra_mask=()):
    """
    Build the recursive 4-copy + blockers backbone.

    Parameters:
      depth: levels of recursion
      cycle_start: which offsets pattern to start with among 4-cycle
      schedule: 'split' or 'after' (blocker arrival vs copies)
      extra_mask: tuple of level indices at which to add a 5th copy

    Returns:
      raw (float) intervals
    """
    T = [(0.0, 1.0)]
    for lvl in range(depth):
        lo = min(l for l, r in T)
        hi = max(r for l, r in T)
        delta = hi - lo

        offs = list(_OFFSETS_CYCLE[(cycle_start + lvl) % 4])
        if lvl in extra_mask:
            # add a small deterministic extra copy at next arithmetic step
            offs.append(max(offs) + 4)

        def _make_copies(src, starts):
            S = []
            for s in starts:
                off = delta * s - lo
                for (l, r) in src:
                    S.append((l + off, r + off))
            return S

        # Use canonical blockers (template 0) by default here; alternates applied by caller
        blockers = [
            (delta * 1, delta * 5),
            (delta * 12, delta * 16),
            (delta * 4, delta * 9),
            (delta * 8, delta * 13),
        ]

        if schedule == 'after':
            T = _make_copies(T, offs) + blockers
        else:  # 'split'
            h = len(offs) // 2
            T = _make_copies(T, offs[:h]) + blockers + _make_copies(T, offs[h:])
    return T

def _apply_blockers_variant(intervals, depth, cycle_start, blockers_idx, schedule, extra_mask):
    """
    Reconstruct a backbone identical to _build_backbone but swapping the blocker template.
    """
    T = [(0.0, 1.0)]
    for lvl in range(depth):
        lo = min(l for l, r in T)
        hi = max(r for l, r in T)
        delta = hi - lo

        offs = list(_OFFSETS_CYCLE[(cycle_start + lvl) % 4])
        if lvl in extra_mask:
            offs.append(max(offs) + 4)

        def _make_copies(src, starts):
            S = []
            for s in starts:
                off = delta * s - lo
                for (l, r) in src:
                    S.append((l + off, r + off))
            return S

        blk_tpl = _BLOCKERS_TEMPLATES[blockers_idx]
        blockers = [(delta * a, delta * b) for (a, b) in blk_tpl]

        if schedule == 'after':
            T = _make_copies(T, offs) + blockers
        else:
            h = len(offs) // 2
            T = _make_copies(T, offs[:h]) + blockers + _make_copies(T, offs[h:])
    return T

# -------------------- Coverage analysis and waves --------------------

def _coverage_cells(intervals):
    """
    Return coverage cells (x_i, x_{i+1}, m_i) where coverage is constant.
    """
    if not intervals:
        return []
    events = []
    for (l, r) in intervals:
        if l < r:
            events.append((l, +1))
            events.append((r, -1))
    events.sort(key=lambda e: (e[0], 0 if e[1] == -1 else 1))
    xs = sorted({x for x, _ in events})
    ev_idx = 0
    cur = 0
    cells = []
    for i in range(len(xs) - 1):
        x = xs[i]
        y = xs[i + 1]
        while ev_idx < len(events) and events[ev_idx][0] == x:
            cur += events[ev_idx][1]
            ev_idx += 1
        if x < y:
            cells.append((x, y, cur))
    return cells

def _try_add_wave_once(current, target_omega, wave_len_set=(2, 3), insert_positions_perc=(1.0, 0.9, 0.75, 0.66, 0.5)):
    """
    Add a single short wave if it improves FirstFit and keeps omega <= target_omega.
    Deterministic selection:
      - choose the longest corridor with coverage <= omega-1
      - try wave lengths in wave_len_set
      - test multiple insertion positions in the arrival sequence
    """
    if not current:
        return None

    base_cols = _firstfit_colors(current)
    base_om = _clique_number(current)
    target_omega = min(target_omega, base_om)

    cells = _coverage_cells(current)
    if not cells:
        return None

    # find longest corridor with coverage <= omega - 1
    L = None
    R = None
    span = -1
    accL = None
    accR = None
    prev_r = None
    for (a, b, m) in cells:
        if m <= base_om - 1:
            if accL is None:
                accL, accR = a, b
            else:
                if prev_r is not None and prev_r == a:
                    accR = b
                else:
                    # finalize
                    if accR - accL > span:
                        L, R, span = accL, accR, accR - accL
                    accL, accR = a, b
        else:
            if accL is not None and accR - accL > span:
                L, R, span = accL, accR, accR - accL
            accL = accR = None
        prev_r = b
    if accL is not None and accR - accL > span:
        L, R, span = accL, accR, accR - accL

    if L is None or R - L <= 0:
        return None

    # build candidate intervals centered in the corridor
    mid = (L + R) // 2
    # candidates' left endpoints around mid
    lefts = sorted(set([L + 1, (2 * L + R) // 3, mid - 1, mid, mid + 1, (L + 2 * R) // 3, max(L + 1, R - 4)]))
    N = len(current)
    probe_positions = sorted(set([max(0, min(N, int(N * p))) for p in insert_positions_perc] + [N]))

    best_seq = None
    best_cols = base_cols
    for wl in wave_len_set:
        for a in lefts:
            b = a + wl
            if a >= b or b > R:
                continue
            cand = (a, b)
            # quick omega check on append
            if _clique_number(current + [cand]) > target_omega:
                continue
            for pos in probe_positions:
                seq = current[:pos] + [cand] + current[pos:]
                if _clique_number(seq) > target_omega:
                    continue
                cols = _firstfit_colors(seq)
                if cols > best_cols:
                    best_cols = cols
                    best_seq = seq
                    # deterministic: return first improving candidate
                    return best_seq
    return None

# -------------------- Frontier-pruning --------------------

def _frontier_prune(intervals, ratio_tol=1e-12, passes=2):
    """
    Deterministic greedy pruning:
      - remove intervals that do not reduce the observed FirstFit/omega ratio
    Priority: try removing shortest intervals last (keep short witnesses).
    """
    cur = list(intervals)
    if not cur:
        return cur
    base_cols = _firstfit_colors(cur)
    base_om = _clique_number(cur) or 1
    base_ratio = base_cols / base_om

    for _ in range(passes):
        changed = True
        while changed:
            changed = False
            order = sorted(range(len(cur)), key=lambda i: (cur[i][1] - cur[i][0], i), reverse=True)
            for idx in order:
                cand = cur[:idx] + cur[idx + 1:]
                if not cand:
                    continue
                om = _clique_number(cand)
                if om == 0:
                    continue
                cols = _firstfit_colors(cand)
                ratio = cols / om
                if ratio + ratio_tol >= base_ratio:
                    cur = cand
                    base_cols, base_om, base_ratio = cols, om, ratio
                    changed = True
                    break
    return cur

# -------------------- Candidate enumeration and selection --------------------

def _evaluate_candidate(depth, cycle_start, schedule, extra_mask, blockers_idx):
    """
    Build candidate backbone with given parameters and compute normalized FF/omega.
    Pattern-aware memoization by signature.
    """
    if not hasattr(_evaluate_candidate, "_cache"):
        _evaluate_candidate._cache = {}
    sig = (depth, cycle_start, schedule, tuple(sorted(extra_mask)), blockers_idx)
    if sig in _evaluate_candidate._cache:
        return _evaluate_candidate._cache[sig]

    raw = _apply_blockers_variant(
        _build_backbone(depth, cycle_start=cycle_start, schedule=schedule, extra_mask=extra_mask),
        depth, cycle_start, blockers_idx, schedule, extra_mask
    )
    norm = _normalize_grid(raw)
    om = _clique_number(norm)
    cols = _firstfit_colors(norm)
    ratio = cols / om if om > 0 else 0.0
    _evaluate_candidate._cache[sig] = (ratio, cols, om, norm)
    return _evaluate_candidate._cache[sig]

# -------------------- Top-level construction --------------------

def construct_intervals():
    """
    Construct a sequence of open intervals maximizing FirstFit/omega.
    Deterministic search over a compact parameter grid + wave augmentation + pruning.
    Returns a list of integer intervals (l, r) on a compact grid.
    """
    # Search parameters (can be tuned)
    depths = [3, 4, 5]
    cycle_starts = [0, 1, 2, 3]                  # rotate A,B,C,D
    schedules = ['split', 'after']
    extra_masks = [(), (0,), (0, 2)]             # no extra, extra at level 0, extras at {0,2}
    blockers_indices = [0, 1]                    # canonical and shifted blockers

    best_ratio = -1.0
    best_cols = -1
    best_om = 1
    best_seq = None

    for d in depths:
        for cs in cycle_starts:
            for sch in schedules:
                for ex in extra_masks:
                    for bi in blockers_indices:
                        ratio, cols, om, seq = _evaluate_candidate(d, cs, sch, ex, bi)
                        # Primary: maximize ratio; tie-break: more colors, then fewer intervals
                        if ratio > best_ratio + 1e-12 or \
                           (abs(ratio - best_ratio) <= 1e-12 and (cols > best_cols or (cols == best_cols and (best_seq is None or len(seq) < len(best_seq))))):
                            best_ratio = ratio
                            best_cols = cols
                            best_om = om
                            best_seq = seq

    if best_seq is None:
        # Fallback: canonical Figure-4 depth=4
        T = [(0.0, 1.0)]
        for _ in range(4):
            lo = min(l for l, r in T)
            hi = max(r for l, r in T)
            delta = hi - lo
            S = []
            for start in (2, 6, 10, 14):
                off = delta * start - lo
                S += [(off + l, off + r) for (l, r) in T]
            S += [
                (delta * 1, delta * 5),
                (delta * 12, delta * 16),
                (delta * 4,  delta * 9),
                (delta * 8,  delta * 13),
            ]
            T = S
        return _normalize_grid(T)

    # Deterministic augmentation with short waves
    current = list(best_seq)
    target_omega = _clique_number(current)
    max_waves = 64
    wave_len_set = (2, 3)
    insert_positions_perc = (1.0, 0.9, 0.75, 0.66, 0.5)

    waves_added = 0
    improved = True
    while improved and waves_added < max_waves:
        improved = False
        candidate = _try_add_wave_once(current, target_omega, wave_len_set=wave_len_set, insert_positions_perc=insert_positions_perc)
        if candidate is not None:
            new_om = _clique_number(candidate)
            new_cols = _firstfit_colors(candidate)
            if new_om <= target_omega and new_cols > _firstfit_colors(current):
                current = candidate
                waves_added += 1
                improved = True

    # Frontier-preserving pruning
    pruned = _frontier_prune(current, ratio_tol=1e-12, passes=2)

    # Normalize again to compact grid and return
    return _normalize_grid(pruned)

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()