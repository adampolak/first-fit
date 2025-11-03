# EVOLVE-BLOCK-START

def construct_intervals(time_budget_per_omega=0.5, max_additions=400):
    """
    Heuristic construction of a sequence of open intervals presented to FirstFit.
    The routine runs a small, time-bounded greedy randomized search for each target
    clique size (omega) and returns the best sequence found (maximizing FirstFit / omega).

    Returns:
      intervals: list of tuples (l, r) representing open intervals in presentation order.
    """
    random.seed(0)

    # Base gadget (Figure 3 from referenced work) - produces FF>OPT on a small instance.
    base = [
        (2.0, 3.0),
        (6.0, 7.0),
        (10.0, 11.0),
        (14.0, 15.0),
        (1.0, 5.0),
        (12.0, 16.0),
        (4.0, 9.0),
        (8.0, 13.0),
    ]

    # Helper: check intersection of open intervals (l1,r1) and (l2,r2)
    def intersects(a, b):
        return (a[0] < b[1]) and (b[0] < a[1])

    # First-Fit simulator (returns colors list and max color)
    def simulate_firstfit(intervals):
        colors = []
        for i, (l, r) in enumerate(intervals):
            used = set()
            # gather colors of previous intervals that intersect this one
            for j in range(i):
                if intersects(intervals[j], (l, r)):
                    used.add(colors[j])
            c = 1
            while c in used:
                c += 1
            colors.append(c)
        return colors, (max(colors) if colors else 0)

    # Clique number (OPT): maximum number of intervals covering a single point
    def compute_clique(intervals):
        events = []
        for (l, r) in intervals:
            # For open intervals, ensure that end events are processed before start events at the same coordinate.
            # Using (x, delta) with delta -1 for end and +1 for start achieves that.
            events.append((l, 1))
            events.append((r, -1))
        events.sort(key=lambda x: (x[0], x[1]))
        cur = 0
        best = 0
        for _, delta in events:
            cur += delta
            if cur > best:
                best = cur
        return best

    # Quick estimate of what color FirstFit would assign to a candidate interval
    # (since FirstFit only uses earlier intervals' colors, we can compute this without full recoloring).
    def candidate_firstfit_color(intervals, colors, cand):
        used = set()
        l, r = cand
        for (iv, c) in zip(intervals, colors):
            if intersects(iv, cand):
                used.add(c)
        c = 1
        while c in used:
            c += 1
        return c

    best_solution = base[:]
    best_colors, best_alg = simulate_firstfit(best_solution)
    best_opt = compute_clique(best_solution)
    best_ratio = (best_alg / best_opt) if best_opt > 0 else 0.0

    # Parameters for candidate generation
    length_choices = [0.3, 0.6, 1.0, 1.6, 2.5, 4.0, 6.0, 10.0]
    jitter_scale = 0.2
    candidate_pool_size = 120  # number of candidate intervals to evaluate per greedy step

    # Try several small target omegas; the algorithm attempts to keep clique <= omega
    for target_w in range(2, 7):  # 2..6
        # start from base gadget scaled a bit to reduce degeneracies
        intervals = [(l * 1.0, r * 1.0) for (l, r) in base]
        colors, alg = simulate_firstfit(intervals)
        opt = compute_clique(intervals)
        if opt > target_w:
            # cannot meet this omega starting from base; skip
            continue

        start_time = time.time()
        steps = 0

        while steps < max_additions and (time.time() - start_time) < time_budget_per_omega:
            steps += 1

            # prepare sample anchor positions: endpoints and midpoints
            endpoints = []
            for (l, r) in intervals:
                endpoints.append(l)
                endpoints.append(r)
            endpoints = sorted(set(endpoints))
            anchors = endpoints[:]
            for i in range(len(endpoints) - 1):
                anchors.append(0.5 * (endpoints[i] + endpoints[i + 1]))
            if not anchors:
                anchors = [0.0]

            cand_list = []
            min_pos = min(endpoints) if endpoints else 0.0
            max_pos = max(endpoints) if endpoints else 1.0
            span = max(1.0, max_pos - min_pos)

            # include a few deterministic 'global' long candidates across the whole span
            global_candidates = [
                (min_pos - 0.2 * span, max_pos + 0.2 * span),
                (min_pos - 0.6 * span, max_pos + 0.2 * span),
            ]

            # generate random and anchored candidates
            for g in global_candidates:
                cand_list.append(g)

            for _ in range(candidate_pool_size):
                base_anchor = random.choice(anchors)
                length = random.choice(length_choices)
                jitter = random.uniform(-jitter_scale, jitter_scale) * max(1.0, length)
                left = base_anchor - length * (0.4 + 0.2 * random.random()) + jitter
                right = left + length
                # Also sometimes create right-anchored candidates
                if random.random() < 0.2:
                    right_anchor = random.choice(anchors)
                    right = right_anchor + random.uniform(-0.1, 0.1) * max(1.0, length)
                    left = right - length
                # Add some random long intervals (bridge many)
                if random.random() < 0.07:
                    left = min_pos - random.uniform(0.1, 0.5) * span
                    right = max_pos + random.uniform(0.1, 0.5) * span
                if left < right:
                    cand_list.append((left, right))

            # Deduplicate and clamp pool
            unique_cands = []
            seen = set()
            for (l, r) in cand_list:
                key = (round(l, 6), round(r, 6))
                if key in seen:
                    continue
                seen.add(key)
                unique_cands.append((l, r))

            # Evaluate candidates: compute predicted color quickly and keep only those that increase max color
            cur_max_color = alg
            promising = []
            for cand in unique_cands:
                # quick neighbor-based color (cheap)
                newcol = candidate_firstfit_color(intervals, colors, cand)
                if newcol > cur_max_color:
                    promising.append((newcol, cand))
            if not promising:
                # no candidate increases FirstFit color under this omega; stop the greedy loop
                break

            # sort by best new color descending; tie-breaker by shorter interval (prefer local gadgets)
            promising.sort(key=lambda x: (-x[0], x[1][1] - x[1][0]))

            accepted = False
            # Try promising candidates in descending order until one is accepted (clique constraint)
            for newcol, cand in promising:
                cand_intervals = intervals + [cand]
                cand_opt = compute_clique(cand_intervals)
                if cand_opt <= target_w:
                    # accept candidate
                    intervals.append(cand)
                    # update colors and alg incrementally: colors remain for earlier intervals; compute new color for appended
                    appended_color = candidate_firstfit_color(intervals[:-1], colors, cand)
                    colors.append(appended_color)
                    alg = max(alg, appended_color)
                    opt = cand_opt
                    accepted = True
                    break
                # otherwise try next promising candidate

            if not accepted:
                # we had promising candidates (would increase FF), but none satisfied the clique constraint.
                # Try to accept a candidate that increases by 1 and keeps clique small by relaxing selection.
                # (This gives a chance for incremental progress.)
                for newcol, cand in sorted(promising, key=lambda x: (x[0], x[1][1] - x[1][0])):
                    cand_intervals = intervals + [cand]
                    cand_opt = compute_clique(cand_intervals)
                    if cand_opt <= target_w:
                        intervals.append(cand)
                        appended_color = candidate_firstfit_color(intervals[:-1], colors, cand)
                        colors.append(appended_color)
                        alg = max(alg, appended_color)
                        opt = cand_opt
                        accepted = True
                        break
                if not accepted:
                    break  # no acceptable candidate this step; stop building for this omega

        # end greedy loop for this omega
        if opt > 0:
            ratio = alg / opt
        else:
            ratio = 0.0

        if ratio > best_ratio:
            best_ratio = ratio
            best_solution = intervals[:]
            best_colors, best_alg = simulate_firstfit(best_solution)
            best_opt = compute_clique(best_solution)

    # final fallback: ensure best_solution is non-empty
    if not best_solution:
        return base[:]
    # Normalize coordinates to tidy integers where possible for readability
    # (scale down very small fractional noise)
    normalized = []
    for (l, r) in best_solution:
        # round small floats to 6 decimals, then if they are near integers convert
        L = round(l, 6)
        R = round(r, 6)
        # promote small fractional parts to a fixed grid to reduce evaluation quirks
        normalized.append((float(L), float(R)))
    return normalized

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()