# EVOLVE-BLOCK-START

def construct_intervals(seed_count=1):
    """
    Adversarial FirstFit generator.
    Input: seed_count (keeps deterministic variations).
    Output: list of (l, r) integer tuples (open intervals) in the order they are presented.
    Guarantees: max number of intervals <= CAP and final clique (opt) <= CLIQUE_LIMIT (10).
    """
    import math
    from bisect import bisect_left, bisect_right
    from collections import defaultdict

    CAP = 9800
    CLIQUE_LIMIT = 10
    rng_seed = int(seed_count) & 0x7FFFFFFF

    # Deterministic RNG (small, fast)
    class SimpleRNG:
        def __init__(self, seed):
            self.x = (seed ^ 0x9e3779b9) & 0xFFFFFFFF
        def randint(self, a, b):
            # xorshift-ish deterministic step
            self.x ^= (self.x << 13) & 0xFFFFFFFF
            self.x ^= (self.x >> 17) & 0xFFFFFFFF
            self.x ^= (self.x << 5) & 0xFFFFFFFF
            return a + (self.x % (b - a + 1))
        def choice(self, seq):
            if not seq:
                return None
            return seq[self.randint(0, len(seq) - 1)]
        def shuffle(self, seq):
            for i in range(len(seq)-1, 0, -1):
                j = self.randint(0, i)
                seq[i], seq[j] = seq[j], seq[i]

    rng = SimpleRNG(rng_seed)

    # Overlap test for open intervals (l1, r1) and (l2, r2)
    def overlaps(a, b):
        return a[0] < b[1] and b[0] < a[1]

    # FirstFit color for a single candidate interval given current intervals with known colors
    def candidate_firstfit_color(intervals, colors, candidate):
        used = set()
        cl, cr = candidate
        for iv, col in zip(intervals, colors):
            if iv[0] < cr and cl < iv[1]:
                used.add(col)
        c = 1
        while c in used:
            c += 1
        return c

    # Sweep-line to compute segment coverage (list of (seg_l, seg_r, count))
    def compute_segments(intervals):
        events = []
        for l, r in intervals:
            # ensure integer endpoints; open intervals: (l, r)
            events.append((l, 1))
            events.append((r, -1))
        if not events:
            return []
        events.sort()
        segments = []
        cur = 0
        last_x = events[0][0]
        i = 0
        n = len(events)
        while i < n:
            x = events[i][0]
            # finalize segment [last_x, x) with count cur
            if x > last_x and cur > 0:
                segments.append((last_x, x, cur))
            # process all events at x
            while i < n and events[i][0] == x:
                cur += events[i][1]
                i += 1
            last_x = x
        # no trailing segment (open intervals)
        return segments

    # Fast max-coverage check when adding one candidate (uses segments recomputed on the fly)
    def max_coverage_with_candidate(intervals, candidate):
        # We build events for intervals + candidate
        events = []
        for l, r in intervals:
            events.append((l, 1))
            events.append((r, -1))
        events.append((candidate[0], 1))
        events.append((candidate[1], -1))
        events.sort()
        cur = 0
        maxc = 0
        i = 0
        n = len(events)
        while i < n:
            x = events[i][0]
            # process all events at x
            while i < n and events[i][0] == x:
                cur += events[i][1]
                i += 1
            if cur > maxc:
                maxc = cur
                # early exit if already exceeding limit
                if maxc > CLIQUE_LIMIT:
                    return maxc
        return maxc

    # Utility: find a "low-coverage" touch point inside interval iv (prefer endpoints + small shift)
    def pick_touch_points_for_interval(iv, segments):
        # prefer right endpoint - 0.5, left endpoint + 0.5, and mid
        l, r = iv
        candidates = []
        if r - l >= 2:
            candidates.append(r - 1)
            candidates.append(l + 1)
        mid = (l + r) // 2
        candidates.append(mid)
        # score them by coverage (we approximate by looking into segments)
        def coverage_at(x):
            for a, b, c in segments:
                if a <= x < b:
                    return c
            return 0
        scored = sorted(candidates, key=lambda x: coverage_at(x))
        return scored

    # Normalize intervals to integer endpoints and r > l
    def normalize_intervals(ints):
        out = []
        for l, r in ints:
            li = int(l)
            ri = int(r)
            if ri <= li:
                ri = li + 1
            out.append((li, ri))
        return out

    # Append helper to add interval and update colors incrementally
    def append_interval(intervals, colors, iv):
        # compute color for iv vs existing intervals
        col = candidate_firstfit_color(intervals, colors, iv)
        intervals.append(iv)
        colors.append(col)
        return col

    # ---------- Build initial sparse anchors ----------
    intervals = []
    colors = []
    # Reserve capacity slices: anchors + adversarial connectors + micro-phase
    # We'll make a moderately large anchor set but not too large to leave room for connectors.
    anchor_count = min(300, max(40, CAP // 40))  # heuristic: ~ CAP/40 or min 40, capped at 300
    gap = 4096  # big spacing so anchors are isolated by default

    base_shift = rng_seed % 37  # small deterministic offset for variety
    for i in range(anchor_count):
        if len(intervals) >= CAP:
            break
        # anchor position widely spaced; tiny intervals of length 1
        pos = (i * gap) + (base_shift % 7) + (i % 3)
        iv = (pos, pos + 1)
        append_interval(intervals, colors, iv)

    # Precompute segments once for selection heuristics
    segments = compute_segments(intervals)

    # Build mapping color -> list of indices (representatives)
    def build_color_reps(colors_list, intervals_list):
        reps = {}
        for idx, col in enumerate(colors_list):
            if col not in reps:
                reps[col] = idx
        return reps

    # ---------- Adversarial connector phase ----------
    # Attempt to force new colors iteratively until we hit capacity or can't find a candidate
    # Limit the number of connector attempts to avoid long run-time
    max_connector_attempts = CAP - len(intervals) - 16  # leave some room for micro-phase
    attempts_done = 0
    stuck_rounds = 0

    while len(intervals) < CAP - 32 and attempts_done < max_connector_attempts and stuck_rounds < 80:
        attempts_done += 1
        if not colors:
            break
        current_max_color = max(colors)
        # small cap: do not try to push colors to absurd numbers beyond a cushion times CLIQUE_LIMIT
        # but allow growth: target up to CLIQUE_LIMIT * 6 or CAP-based
        reasonable_upper = max(CLIQUE_LIMIT * 6, 60)
        if current_max_color >= reasonable_upper:
            break

        reps = build_color_reps(colors, intervals)
        # We must pick for each color 1..current_max_color some representative interval to intersect
        # If some color is missing (shouldn't happen), skip;
        if any(c not in reps for c in range(1, current_max_color + 1)):
            # rebuild mapping and continue
            reps = build_color_reps(colors, intervals)
            if any(c not in reps for c in range(1, current_max_color + 1)):
                break

        # collect representative intervals
        rep_indices = [reps[c] for c in range(1, current_max_color + 1)]
        rep_intervals = [intervals[idx] for idx in rep_indices]

        # compute segments for current state once
        segments = compute_segments(intervals)

        # For better chances avoid representatives clustered around same spot:
        # reorder reps by the segment coverage at their midpoints (prefer low coverage first)
        def coverage_at_point(x):
            for a, b, c in segments:
                if a <= x < b:
                    return c
            return 0
        rep_points = []
        for iv in rep_intervals:
            l, r = iv
            mid = (l + r) // 2
            rep_points.append((coverage_at_point(mid), iv))
        rep_points.sort(key=lambda x: x[0])  # low coverage first
        ordered_reps = [iv for _, iv in rep_points]

        # Try a few variants to find a feasible connector forcing a new color
        found = False
        # We'll try simple variants: touching endpoints (right-1/left+1), midpoints, and sparse subset spanning
        variants = []
        # full-span variant (touch near endpoints of each rep)
        variants.append(("full_touch", ordered_reps))
        # sparse variant: sample every ceil(step) rep to reduce overlap density
        step = max(1, len(ordered_reps) // 8)
        variants.append(("sparse", ordered_reps[::step]))
        # reversed order
        variants.append(("rev", list(reversed(ordered_reps))))
        # chunked spans (try multiple windows)
        chunk = max(1, len(ordered_reps) // 4)
        for cstart in range(0, len(ordered_reps), chunk):
            variants.append(("chunk", ordered_reps[cstart:cstart + chunk]))

        # Try each variant with some touch-point heuristics
        for vname, reps_try in variants:
            if found:
                break
            if not reps_try:
                continue
            # choose candidate touch points for each representative: try right-edge, left-edge, mid variants
            # Build arrays of touch candidates per rep (1..3 choices)
            touch_choices = []
            for iv in reps_try:
                l, r = iv
                choices = []
                # prefer to intersect near an endpoint to reduce length
                if r - l >= 2:
                    choices.append(r - 1)
                    choices.append(l + 1)
                choices.append((l + r) // 2)
                # optionally include leftmost or rightmost extremes
                choices = list(dict.fromkeys(choices))  # unique while preserving order
                touch_choices.append(choices)

            # We'll try a few combinations greedily: pick best point per rep by low coverage, then fallback to some random mixes
            # 1) deterministic best-per-rep
            sel_points = []
            for choices in touch_choices:
                # score by current coverage
                best = None
                best_cov = 10**9
                for p in choices:
                    cov = 0
                    for a, b, c in segments:
                        if a <= p < b:
                            cov = c
                            break
                    if cov < best_cov:
                        best_cov = cov
                        best = p
                sel_points.append(best)
            if sel_points:
                s = min(sel_points) - 1
                e = max(sel_points) + 1
                if e <= s:
                    e = s + 1
                candidate = (s, e)
                # quick check: candidate must intersect at least one interval of each color in reps_try
                ok_intersect_all = True
                for iv in reps_try:
                    if not overlaps(iv, candidate):
                        ok_intersect_all = False
                        break
                if not ok_intersect_all:
                    # try to extend slighty to ensure intersection
                    candidate = (s - 1, e + 1)
                # compute what color FirstFit would assign
                target_color = candidate_firstfit_color(intervals, colors, candidate)
                if target_color == current_max_color + 1:
                    # ensure max coverage remains <= CLIQUE_LIMIT
                    maxcov = max_coverage_with_candidate(intervals, candidate)
                    if maxcov <= CLIQUE_LIMIT:
                        # accept
                        append_interval(intervals, colors, candidate)
                        found = True
                        stuck_rounds = 0
                        break

            # 2) random mixing tries (deterministic RNG)
            mix_tries = 6
            for _ in range(mix_tries):
                sel_points = []
                for choices in touch_choices:
                    sel_points.append(rng.choice(choices))
                s = min(sel_points) - rng.randint(0, 2) - 1
                e = max(sel_points) + rng.randint(0, 2) + 1
                if e <= s:
                    e = s + 1
                candidate = (s, e)
                # check overlap with reps_try
                ok_intersect_all = all(overlaps(iv, candidate) for iv in reps_try)
                if not ok_intersect_all:
                    continue
                target_color = candidate_firstfit_color(intervals, colors, candidate)
                if target_color == current_max_color + 1:
                    maxcov = max_coverage_with_candidate(intervals, candidate)
                    if maxcov <= CLIQUE_LIMIT:
                        append_interval(intervals, colors, candidate)
                        found = True
                        stuck_rounds = 0
                        break

        if not found:
            stuck_rounds += 1
        # update segments occasionally for heuristics
        if len(intervals) % 50 == 0:
            segments = compute_segments(intervals)

    # ---------- If adversary phase stalls, fill with safe micro patterns ----------
    # We'll add fractional-window micro-blocks and towers that are designed to avoid increasing clique.
    def micro_phase_fill(intervals, colors, cap_left):
        if cap_left <= 0:
            return
        lo = 0
        hi = 1
        if intervals:
            lo = min(l for l, _ in intervals)
            hi = max(r for _, r in intervals)
        span = max(1, hi - lo)
        # four window fractions (offsets chosen to diversify)
        windows = [(0.06, 0.14), (0.26, 0.34), (0.48, 0.62), (0.74, 0.86)]
        # sample a small seed from existing intervals (or create synthetic if none)
        seed = intervals[:min(32, max(8, len(intervals) // 200))] if intervals else [(0, 1)]
        if not seed:
            seed = [(0, 1)]
        # create blocks scaled to windows
        blocks = []
        for (fa, fb) in windows:
            win_lo = lo + int(round(fa * span))
            base_shift = win_lo - min(l for l, _ in seed)
            # translate seed
            block = [(l + base_shift, r + base_shift) for (l, r) in seed]
            blocks.append(block)
        # interleave blocks
        micro = []
        maxlen = max(len(b) for b in blocks)
        for i in range(maxlen):
            for blk in blocks:
                if i < len(blk):
                    micro.append(blk[i])
                    if len(micro) >= cap_left:
                        break
            if len(micro) >= cap_left:
                break
        # add a few connectors across windows (safe)
        connectors = []
        connectors.append((lo + int(0.02 * span), lo + int(0.18 * span)))
        connectors.append((lo + int(0.36 * span), lo + int(0.62 * span)))
        connectors.append((lo + int(0.68 * span), lo + int(0.94 * span)))
        # append micro and connectors, validating clique bound
        for iv in micro + connectors:
            if len(intervals) >= CAP:
                break
            if max_coverage_with_candidate(intervals, iv) <= CLIQUE_LIMIT:
                append_interval(intervals, colors, iv)
        return

    room = CAP - len(intervals)
    # run a couple of micro fills to top up safely
    micro_rounds = 6
    for _ in range(micro_rounds):
        room = CAP - len(intervals)
        if room <= 0:
            break
        micro_phase_fill(intervals, colors, min(512, room))

    # final small tower-and-cap fill (short intervals likely to touch many active colors)
    def tower_fill(intervals, colors):
        room = CAP - len(intervals)
        if room <= 0:
            return
        lo = min((l for l, _ in intervals), default=0)
        hi = max((r for _, r in intervals), default=1)
        span = max(1, hi - lo)
        towers = [(0.10, 0.18), (0.29, 0.37), (0.46, 0.54), (0.63, 0.71), (0.80, 0.88)]
        seed = intervals[:min(32, max(8, len(intervals) // 200))] if intervals else [(0,1)]
        for fa, fb in towers:
            if len(intervals) >= CAP:
                break
            win_lo = lo + int(round(fa * span))
            base_shift = win_lo - min(l for l, _ in seed)
            block = []
            eps = max(1, span // 512)
            for idx, (l, r) in enumerate(seed):
                mid = (l + r) // 2
                L = mid + base_shift + (idx % 3)
                R = L + eps
                block.append((L, R))
            # interleave block into intervals as long as safe
            for iv in block:
                if len(intervals) >= CAP:
                    break
                if max_coverage_with_candidate(intervals, iv) <= CLIQUE_LIMIT:
                    append_interval(intervals, colors, iv)

    tower_fill(intervals, colors)

    # Final capacity guard and normalization
    if len(intervals) > CAP:
        intervals = intervals[:CAP]
        colors = colors[:CAP]

    intervals = normalize_intervals(intervals)

    # ensure r > l and integer endpoints
    if len(intervals) > CAP:
        intervals = intervals[:CAP]

    return intervals

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()