# EVOLVE-BLOCK-START

from functools import lru_cache

# ------------------ low-level utilities ------------------

def _normalize_grid(intervals):
    """
    Normalize endpoints to a compact integer grid while preserving order (open intervals).
    Each unique endpoint is mapped to an increasing even integer.
    """
    if not intervals:
        return []
    pts = sorted(set(x for seg in intervals for x in seg))
    coord = {e: 2 * i for i, e in enumerate(pts)}
    return [(coord[l], coord[r]) for (l, r) in intervals]


def _signature(intervals):
    """
    Canonical pattern signature used for memoization:
    - normalize endpoints to ranks (even integers)
    - return tuple of ints for hashing
    """
    Tn = _normalize_grid(intervals)
    flat = []
    for (l, r) in Tn:
        flat.append(l)
        flat.append(r)
    return tuple(flat)


def _ff_detail(intervals):
    """
    FirstFit simulation for open intervals in arrival order.
    Returns:
      colors_used, color_assignment, new_color_indices (indices at which a new color opened)
    Optimization: track last_end per color (sufficient for open intervals with arrival order).
    """
    last_end = []
    color_of = []
    new_color_indices = []
    for idx, (l, r) in enumerate(intervals):
        placed = False
        for c in range(len(last_end)):
            if l >= last_end[c]:
                last_end[c] = r
                color_of.append(c)
                placed = True
                break
        if not placed:
            last_end.append(r)
            color_of.append(len(last_end) - 1)
            new_color_indices.append(idx)
    return len(last_end), color_of, new_color_indices


def _omega_open(intervals):
    """
    Compute clique number (omega) for open intervals via sweep:
    process right endpoints before left endpoints at the same coordinate.
    """
    events = []
    for (l, r) in intervals:
        if l < r:
            events.append((l, +1))
            events.append((r, -1))
    # right(-1) before left(+1)
    events.sort(key=lambda e: (e[0], 0 if e[1] == -1 else 1))
    cur = best = 0
    for _, t in events:
        cur += t
        if cur > best:
            best = cur
    return best


# Cache evaluation by signature to accelerate exploration
_eval_cache = {}

def _evaluate(intervals):
    """
    Evaluate instance: return (ratio, omega, colors, n, color_of, pivots).
    Memoized by signature.
    """
    if not intervals:
        return (-1.0, 0, 0, 0, [], [])
    sig = _signature(intervals)
    if sig in _eval_cache:
        return _eval_cache[sig]
    om = _omega_open(intervals)
    if om <= 0:
        res = (-1.0, 0, 0, len(intervals), [], [])
        _eval_cache[sig] = res
        return res
    cols, color_of, pivots = _ff_detail(intervals)
    ratio = cols / om
    res = (ratio, om, cols, len(intervals), color_of, pivots)
    _eval_cache[sig] = res
    return res


# ------------------ construction primitives ------------------

def _make_copies(T, offsets, delta, lo, center, translation, reverse_alt=False, interleave='block'):
    """
    Create copy-lists of T at offsets and compose according to interleave strategy.
    """
    copy_lists = []
    for idx, start in enumerate(offsets):
        if translation == 'left':
            off = delta * start - lo
        else:
            off = delta * start - center
        seq = T if not (reverse_alt and (idx % 2 == 1)) else list(reversed(T))
        copy_lists.append([(l + off, r + off) for (l, r) in seq])

    # Compose copies
    if interleave == 'block':
        S_copies = []
        for lst in copy_lists:
            S_copies.extend(lst)
    else:  # 'zip'
        S_copies = []
        if copy_lists:
            m = len(copy_lists[0])
            for j in range(m):
                for lst in copy_lists:
                    S_copies.append(lst[j])
    return S_copies


def _scale_blockers(blockers, delta, center, anchor='left'):
    """
    Scale blocker template by delta. Anchor 'left' uses raw multipliers; 'center' translates by center.
    """
    S = []
    for (a, b) in blockers:
        if anchor == 'left':
            S.append((delta * a, delta * b))
        else:
            S.append((delta * a - center, delta * b - center))
    return S


def _add_waves(delta, windows, spec):
    """
    Deterministically generate wave intervals within given windows.
    spec is a dict:
      kind: 'short' | 'long' | 'both'
      short_len: length as fraction of delta (default 0.03)
      long_len: length as fraction of delta (default 0.08)
      density: number of intervals per window (default 6)
      gap: fractional gap step between starts (default 0.12)
    Returns list of intervals.
    """
    if not windows or not spec or spec.get('kind') == 'none':
        return []

    kind = spec.get('kind', 'short')
    short_len = float(spec.get('short_len', 0.03)) * delta
    long_len = float(spec.get('long_len', 0.08)) * delta
    density = int(spec.get('density', 6))
    gap = float(spec.get('gap', 0.12)) * delta

    waves = []
    length_choices = []
    if kind in ('short', 'both'):
        length_choices.append(short_len)
    if kind in ('long', 'both'):
        length_choices.append(long_len)
    if not length_choices:
        return []

    for w_idx, (L, R) in enumerate(windows):
        span = max(0.0, R - L)
        if span <= 0:
            continue
        # start positions staggered deterministically
        base = L + 0.02 * span
        step = max(gap, span / max(1, density + 2))
        for i in range(density):
            length = length_choices[(w_idx + i) % len(length_choices)]
            s = base + i * step
            if s + length < R:
                waves.append((s, s + length))
    return waves


def _build_level(T, starts, blockers, translation, blocker_anchor,
                 schedule, interleave, reverse_alt, wavespec, corridor_windows):
    """
    Build one expansion level from T using copies + blockers + optional waves.
    """
    lo = min(l for l, r in T)
    hi = max(r for l, r in T)
    delta = hi - lo
    center = (lo + hi) / 2.0

    S_copies = _make_copies(T, starts, delta, lo, center, translation,
                            reverse_alt=reverse_alt, interleave=interleave)
    S_blockers = _scale_blockers(blockers, delta, center, anchor=blocker_anchor)
    S_waves = _add_waves(delta, corridor_windows, wavespec)

    if schedule == 'before':
        S = S_blockers + S_copies + S_waves
    elif schedule == 'split' and interleave == 'block':
        h = max(1, len(starts) // 2)
        # reconstruct halves from S_copies by block grouping
        grouped = []
        base = 0
        size = len(T)
        for _ in range(len(starts)):
            grouped.append(S_copies[base:base + size])
            base += size
        first_half = []
        for i in range(h):
            first_half.extend(grouped[i])
        second_half = []
        for i in range(h, len(starts)):
            second_half.extend(grouped[i])
        S = first_half + S_blockers + S_waves + second_half
    elif schedule == 'mix':
        # interleave copies and blockers round-robin, then waves
        S = []
        i = j = 0
        while i < len(S_copies) or j < len(S_blockers):
            if i < len(S_copies):
                S.append(S_copies[i]); i += 1
            if j < len(S_blockers):
                S.append(S_blockers[j]); j += 1
        S += S_waves
    else:
        # 'after'
        S = S_copies + S_blockers + S_waves
    return S


def _build_pattern(k, base_seed, starts, blockers, translation='left',
                   blocker_anchor='left', schedule='after', interleave='block',
                   reverse_alt=False, wavespec=None, corridor='standard', extra_caps=False):
    """
    Expand base_seed k times using the copy + blocker + wave scheme.
    corridor: choose among predefined corridor windows scaled by delta.
    extra_caps: append extra long "cap" intervals per level to couple towers.
    """
    T = list(base_seed)
    for _ in range(k):
        # choose corridor windows in multipliers of delta
        if corridor == 'tight':
            corridor_windows = [(4.5, 7.5), (9.0, 11.0)]
        elif corridor == 'wide':
            corridor_windows = [(3.5, 8.5), (8.5, 13.5)]
        else:
            corridor_windows = [(4.0, 7.0), (8.0, 12.0)]
        # scale corridors by delta later in _add_waves using (L*delta, R*delta)
        # But _add_waves expects absolute coordinates; we need delta at this level.
        # We'll pass absolute windows here:
        lo = min(l for l, r in T)
        hi = max(r for l, r in T)
        delta = hi - lo
        abs_windows = [(delta * a, delta * b) for (a, b) in corridor_windows]
        T = _build_level(
            T, starts, blockers, translation, blocker_anchor,
            schedule, interleave, reverse_alt, wavespec, abs_windows
        )
        if extra_caps:
            # Two additional caps per level across wide spans (do not depend on anchor)
            lo2 = min(l for l, r in T)
            hi2 = max(r for l, r in T)
            d2 = hi2 - lo2
            # Caps sized and placed to touch many towers but avoid expanding omega too much
            caps = [(d2 * 2.0, d2 * 10.0), (d2 * 6.0, d2 * 14.0)]
            T.extend(caps)
    return T


# ------------------ pruning ------------------

def _two_stage_prune(intervals, base_ratio, pin_pivots=True, max_rounds=2):
    """
    Two-stage deterministic pruning.
    Stage 1: preserve frontier witness (pivots that opened new colors).
    Stage 2: aggressive pruning of any non-pinned intervals, accepting removal if ratio >= base_ratio.
    """
    if not intervals:
        return intervals

    # Establish base pivots
    _, _, _, _, _, pivots = _evaluate(intervals)
    pinned = set(pivots) if pin_pivots else set()

    def try_prune(T, allow_pinned=False):
        cur = list(T)
        # Order: prefer removing longer intervals first (and later arrivals earlier)
        lens = [(i, cur[i][1] - cur[i][0]) for i in range(len(cur))]
        order = [i for i, _ in sorted(lens, key=lambda x: (-x[1], -x[0]))]
        changed = True
        while changed:
            changed = False
            for idx in order:
                if idx >= len(cur):
                    continue
                if (not allow_pinned) and idx in pinned:
                    continue
                cand = cur[:idx] + cur[idx+1:]
                ratio, om, cols, _, _, _ = _evaluate(cand)
                if om > 0 and ratio >= base_ratio:
                    cur = cand
                    # recompute order positions if structural change
                    lens = [(i, cur[i][1] - cur[i][0]) for i in range(len(cur))]
                    order = [i for i, _ in sorted(lens, key=lambda x: (-x[1], -x[0]))]
                    changed = True
                    break
        return cur

    # Stage 1: protected pivots
    T1 = try_prune(intervals, allow_pinned=False)
    # Stage 2: a couple of aggressive passes with pinned protection
    T2 = T1
    for _ in range(max_rounds):
        T2 = try_prune(T2, allow_pinned=False)
    return T2


# ------------------ top-level construction ------------------

def construct_intervals(iterations=4):
    """
    Build a sequence of open intervals aiming to maximize FirstFit/omega.
    Deterministic sweep over:
      - depths k ∈ {3,4,5} (plus 'iterations' if in set),
      - offset cycles A/B/C/D,
      - blocker templates (canonical + variants + caps),
      - schedules and interleave modes,
      - translation/anchor modes,
      - corridor windows and wave types,
      - extra caps on/off.
    Selects best ratio and then applies two-stage pruning. Falls back to canonical witness.
    Returns normalized integer endpoints.
    """

    # Tiling patterns (A/B/C/D)
    offset_cycles = [
        (2, 6, 10, 14),
        (1, 5, 9, 13),
        (3, 7, 11, 15),
        (0, 4, 8, 12),
    ]

    # Blocker templates: canonical + variants + cap-heavy
    blocker_templates = [
        # Canonical four-blocker (Figure-4)
        ((1, 5), (12, 16), (4, 9), (8, 13)),
        # Slightly shifted mid couplers
        ((1, 5), (11, 15), (4, 9), (8, 13)),
        # Tighter middle
        ((2, 6), (12, 16), (4, 8.5), (9.5, 13)),
        # Asymmetric pushers
        ((1, 6), (11, 16), (3, 7), (9, 13)),
        # Cap-augmented (interpreted as 6 blockers)
        ((1, 5), (12, 16), (4, 9), (8, 13), (2, 14), (6, 10)),
    ]

    # Base seeds (keep clique modest)
    base_seeds = [
        [(0.0, 1.0)],
        [(0.0, 1.0), (2.0, 3.0)],
    ]

    depths = sorted({d for d in (3, 4, 5, iterations) if d in (3, 4, 5)})
    if not depths:
        depths = [4]

    schedules = ['after', 'before', 'split', 'mix']
    interleaves = ['block', 'zip']
    translations = ['left', 'center']
    anchors = ['left', 'center']
    corridors = ['standard', 'tight', 'wide']
    wave_specs = [
        {'kind': 'none'},
        {'kind': 'short', 'short_len': 0.03, 'density': 6, 'gap': 0.10},
        {'kind': 'long',  'long_len': 0.08, 'density': 4, 'gap': 0.14},
        {'kind': 'both',  'short_len': 0.03, 'long_len': 0.08, 'density': 6, 'gap': 0.12},
    ]
    cap_flags = [False, True]
    reverse_flags = [False, True]

    # Size guard
    MAX_INTERVALS = 2200

    best = None  # (ratio, om, cols, n, intervals)

    for base in base_seeds:
        for k in depths:
            # rough growth check: (~copies^k)*|base| + blockers per level
            for starts in offset_cycles:
                for blockers in blocker_templates:
                    # estimate copies and blockers count
                    copies = len(starts)
                    blk_count = len(blockers)
                    est = (copies ** k) * len(base) + blk_count * k
                    if est > MAX_INTERVALS:
                        continue
                    for schedule in schedules:
                        for inter in interleaves:
                            for tr in translations:
                                for anch in anchors:
                                    for corr in corridors:
                                        for wavespec in wave_specs:
                                            for caps in cap_flags:
                                                for rev in reverse_flags:
                                                    T = _build_pattern(
                                                        k=k, base_seed=base, starts=starts, blockers=blockers,
                                                        translation=tr, blocker_anchor=anch, schedule=schedule,
                                                        interleave=inter, reverse_alt=rev, wavespec=wavespec,
                                                        corridor=corr, extra_caps=caps
                                                    )
                                                    ratio, om, cols, n, _, _ = _evaluate(T)
                                                    # Skip degenerate
                                                    if om == 0 or n == 0:
                                                        continue
                                                    # Select best ratio, then smaller n, then more colors
                                                    cand = (ratio, om, cols, n, T)
                                                    if best is None:
                                                        best = cand
                                                    else:
                                                        if cand[0] > best[0] + 1e-9:
                                                            best = cand
                                                        elif abs(cand[0] - best[0]) <= 1e-9:
                                                            if cand[3] < best[3] or (cand[3] == best[3] and cand[2] > best[2]):
                                                                best = cand

    # Fallback: canonical Figure-4 (depth 4), known FF≈13, ω≈5
    if best is None or best[0] <= 0:
        T = [(0.0, 1.0)]
        for _ in range(4):
            lo = min(l for l, r in T)
            hi = max(r for l, r in T)
            delta = hi - lo
            S = []
            for start in (2, 6, 10, 14):
                off = delta * start - lo
                S += [(l + off, r + off) for (l, r) in T]
            S += [(delta * 1, delta * 5), (delta * 12, delta * 16), (delta * 4, delta * 9), (delta * 8, delta * 13)]
            T = S
        return _normalize_grid(T)

    # Two-stage pruning around the best ratio
    base_ratio = best[0]
    T_best = best[4]
    T_pruned = _two_stage_prune(T_best, base_ratio=base_ratio, pin_pivots=True, max_rounds=2)

    # Ensure final ratio not worse than base (otherwise use unpruned)
    r_p, om_p, cols_p, _, _, _ = _evaluate(T_pruned)
    if om_p > 0 and r_p + 1e-12 >= base_ratio:
        return _normalize_grid(T_pruned)
    else:
        return _normalize_grid(T_best)

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()