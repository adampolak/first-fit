# EVOLVE-BLOCK-START

from math import gcd

# ---------- Core utilities: geometry, FirstFit, clique ----------

def overlaps(a, b):
    (l1, r1), (l2, r2) = a, b
    return max(l1, l2) < min(r1, r2)

def firstfit_colors(intervals):
    """
    Simulate FirstFit coloring for the given arrival order.
    Uses per-color conflict checks (robust to arbitrary arrival order).
    """
    colors = []
    for iv in intervals:
        placed = False
        for c in colors:
            # check conflict against all intervals in color c
            ok = True
            for u in c:
                if overlaps(iv, u):
                    ok = False
                    break
            if ok:
                c.append(iv)
                placed = True
                break
        if not placed:
            colors.append([iv])
    return len(colors)

def clique_number(intervals):
    """
    Compute omega = max number of open intervals covering a single point.
    Sweep with right endpoints processed before left endpoints at ties.
    """
    ev = []
    for (l, r) in intervals:
        if l < r:
            ev.append((l, +1))
            ev.append((r, -1))
    ev.sort(key=lambda x: (x[0], 0 if x[1] == -1 else 1))
    cur = best = 0
    for _, t in ev:
        cur += t
        if cur > best:
            best = cur
    return best

def normalize_to_grid(intervals):
    """
    Map unique endpoints to a compact even-integer grid while preserving order.
    """
    if not intervals:
        return []
    pts = sorted({x for seg in intervals for x in seg})
    mp = {}
    cur = 0
    for x in pts:
        mp[x] = cur
        cur += 2
    return [(mp[l], mp[r]) for (l, r) in intervals]

def normalize_and_reduce_gcd(intervals):
    """
    Normalize to even-integer grid and divide by global gcd to keep footprint small.
    """
    if not intervals:
        return []
    L = normalize_to_grid(intervals)
    g = 0
    for a, b in L:
        g = gcd(g, abs(a))
        g = gcd(g, abs(b))
    if g > 1:
        L = [(a // g, b // g) for (a, b) in L]
    return L

# ---------- Blueprint: 4-copy expansion with blockers ----------

def expand_once(seed, offsets, blockers):
    """
    One expansion level:
      - Place |offsets| translated copies of seed by multipliers in 'offsets'.
      - Add long 'blockers' given as multipliers (a, b) of delta.
    """
    if not seed:
        return []
    lo = min(l for l, r in seed)
    hi = max(r for l, r in seed)
    delta = hi - lo
    S = []
    for s in offsets:
        off = delta * s - lo
        for (l, r) in seed:
            S.append((l + off, r + off))
    for (a, b) in blockers:
        S.append((delta * a, delta * b))
    return S

def build_pattern(depth, offsets, blockers, schedule='after', add_extra_first=False, normalize_schedule='final'):
    """
    Build a multi-level pattern with:
      - depth in {3,4,5} explored by search
      - offsets: tuple of four (or more) start multipliers
      - blockers: list of four connector multipliers
      - schedule: 'before'|'after'|'split' (blockers placement within level)
      - add_extra_first: add a fifth copy at level 0 (offsets + {last+4})
      - normalize_schedule: 'none'|'each'|'every2'|'final'
    """
    T = [(0.0, 1.0)]
    for lvl in range(depth):
        offs = list(offsets)
        if add_extra_first and lvl == 0:
            offs = offs + [offs[-1] + 4]
        # choose schedule
        lo = min(l for l, r in T)
        hi = max(r for l, r in T)
        delta = hi - lo

        def copies_list(sub_offsets):
            S = []
            base_lo = lo
            for s in sub_offsets:
                off = delta * s - base_lo
                for (l, r) in T:
                    S.append((l + off, r + off))
            return S

        block_list = [(delta * a, delta * b) for (a, b) in blockers]

        if schedule == 'before':
            S = list(block_list) + copies_list(offs)
        elif schedule == 'after':
            S = copies_list(offs) + list(block_list)
        else:  # 'split'
            h = len(offs) // 2
            S = copies_list(offs[:h]) + list(block_list) + copies_list(offs[h:])

        T = S

        # optional mid-level normalization to perturb tie structure deterministically
        if normalize_schedule == 'each':
            T = normalize_to_grid(T)
        elif normalize_schedule == 'every2' and (lvl % 2 == 1):
            T = normalize_to_grid(T)

    # final normalization policy
    if normalize_schedule in ('final', 'each', 'every2'):
        T = normalize_and_reduce_gcd(T)
    return T

# ---------- Deterministic frontier-preserving pruning ----------

def prune_preserve_colors_and_omega(intervals):
    """
    Remove intervals if and only if both FirstFit colors and omega remain identical.
    Greedy deterministic pass, removing longer intervals first.
    """
    cur = list(intervals)
    if not cur:
        return cur
    base_cols = firstfit_colors(cur)
    base_om = clique_number(cur)
    def length(iv): return iv[1] - iv[0]
    changed = True
    while changed:
        changed = False
        order = sorted(range(len(cur)), key=lambda i: (-length(cur[i]), i))
        for i in order:
            cand = cur[:i] + cur[i+1:]
            if not cand:
                continue
            c = firstfit_colors(cand)
            o = clique_number(cand)
            if c == base_cols and o == base_om:
                cur = cand
                changed = True
                break
    return cur

# ---------- Search space enumerator and selector ----------

def construct_intervals():
    """
    Enumerate a deterministic family of blueprints (depth/offsets/blockers/schedule/normalization),
    select the one that maximizes FirstFit/omega, tie-breaking by fewer intervals and then
    by larger FirstFit colors. Finally prune deterministically while preserving the frontier.
    """
    # Offsets sets (four-copy placements)
    offset_sets = [
        (2, 6, 10, 14),   # A
        (1, 5, 9, 13),    # B
        (3, 7, 11, 15),   # C
        (0, 4, 8, 12),    # D
    ]

    # Blocker templates (four connectors)
    blocker_sets = [
        ((1, 5), (12, 16), (4, 9), (8, 13)),    # Template A (canonical)
        ((0, 4), (11, 15), (3, 8), (7, 12)),    # Template B
        ((1, 5), (9, 13), (5, 9), (12, 16)),    # Template C (shifted pairs)
    ]

    depths = [3, 4, 5]
    schedules = ['after', 'before', 'split']
    normalize_policies = ['final', 'each', 'every2']  # 'none' tends to grow footprint; keep compact
    extra_first_choices = [False, True]

    best_T = None
    best_ratio = -1.0
    best_n = None
    best_cols = -1
    best_om = 1

    for k in depths:
        for offs in offset_sets:
            for blks in blocker_sets:
                for sch in schedules:
                    for normp in normalize_policies:
                        for ex in extra_first_choices:
                            T = build_pattern(
                                depth=k,
                                offsets=offs,
                                blockers=blks,
                                schedule=sch,
                                add_extra_first=ex,
                                normalize_schedule=normp
                            )
                            if not T:
                                continue
                            om = clique_number(T)
                            if om == 0:
                                continue
                            cols = firstfit_colors(T)
                            ratio = cols / om
                            n = len(T)
                            better = False
                            if ratio > best_ratio + 1e-12:
                                better = True
                            elif abs(ratio - best_ratio) <= 1e-12:
                                if best_n is None or n < best_n - 0:
                                    better = True
                                elif n == best_n and cols > best_cols:
                                    better = True
                            if better:
                                best_ratio = ratio
                                best_T = T
                                best_n = n
                                best_cols = cols
                                best_om = om

    # Fallback to canonical if search fails (should not happen)
    if best_T is None:
        best_T = build_pattern(
            depth=4,
            offsets=(2, 6, 10, 14),
            blockers=((1, 5), (12, 16), (4, 9), (8, 13)),
            schedule='after',
            add_extra_first=False,
            normalize_schedule='final'
        )

    # Deterministic pruning preserving both FF colors and omega
    pruned = prune_preserve_colors_and_omega(best_T)
    # Ensure normalization to compact grid and return
    return normalize_and_reduce_gcd(pruned)

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()