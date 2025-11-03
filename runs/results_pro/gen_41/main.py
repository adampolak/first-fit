# EVOLVE-BLOCK-START

def construct_intervals():
    """
    Build a recursive 4-layer interval gadget with diversified translations/bridges,
    then adversarially order intervals via a greedy mex scheduler to inflate
    FirstFit color usage while keeping the clique small.

    Returns:
      List of open intervals (l, r) presented to FirstFit in the chosen order.
    """
    import math

    # ------------------------ Basic geometry utilities ------------------------
    def overlap(a, b):
        # Open-interval overlap
        return (a[0] < b[1]) and (b[0] < a[1])

    def clique_size(intervals):
        # Maximum number of intervals covering a single point (order-independent)
        events = []
        for (l, r) in intervals:
            events.append((l, 1))   # start
            events.append((r, -1))  # end
        # Process ends before starts at the same coordinate for open intervals
        events.sort(key=lambda x: (x[0], x[1]))
        cur = 0
        best = 0
        best_x = None
        for x, d in events:
            cur += d
            if cur > best:
                best = cur
                best_x = x
        return best, (best_x if best_x is not None else 0.0)

    # ------------------------ Recursive backbone builder ----------------------
    def build_backbone(depth=4):
        """
        Four-level recursive expansion, cycling start patterns and bridge-sets
        with a gentle gamma-schedule to diversify but preserve small omega.
        """
        # Multi-pattern (translation starts) cycle
        START_PATTERNS = [
            (2, 6, 10, 14),
            (3, 7, 11, 15),
            (4, 8, 12, 16),
            (5, 9, 13, 17),
        ]
        # Bridge-set cycle (scaled by current span per level)
        BRIDGE_SETS = [
            [(1, 5), (12, 16), (4, 9), (8, 13)],
            [(2, 6), (11, 15), (5, 10), (9, 14)],
            [(3, 7), (10, 14), (6, 11), (7, 12)],
            [(4, 8), (13, 17), (7, 12), (8, 13)],
        ]
        # Gentle gamma-schedule per level
        GAMMA = [1.00, 1.20, 0.95, 1.15]

        # Two short disjoint seeds produce rich structure early, but keep omega modest
        T = [(0.0, 1.0), (2.0, 3.0)]

        for i in range(depth):
            lo = min(l for l, r in T)
            hi = max(r for l, r in T)
            base_span = (hi - lo)
            span = base_span * GAMMA[i % len(GAMMA)]

            starts = START_PATTERNS[i % len(START_PATTERNS)]
            bridges = BRIDGE_SETS[i % len(BRIDGE_SETS)]
            S = []
            # Replicate T using the chosen start offsets
            for s in starts:
                offset = span * s - lo
                for (l, r) in T:
                    S.append((l + offset, r + offset))
            # Inject bridges
            for (a, b) in bridges:
                S.append((span * a, span * b))
            T = S

        return T

    def normalize(intervals, target_span=2000.0):
        """
        Shift to positive, rescale so total span ~ target_span, and keep
        fractional coordinates to avoid accidental boundary coincidences.
        """
        lo = min(l for l, r in intervals)
        hi = max(r for l, r in intervals)
        span = hi - lo if hi > lo else 1.0
        # Shift so min starts above 1.0
        shift = (-lo + 1.0) if lo <= 0.0 else 0.0
        scale = target_span / span
        out = []
        for (l, r) in intervals:
            L = (l + shift) * scale
            R = (r + shift) * scale
            # Keep as floats; ensure open intervals remain valid
            if R <= L:
                R = L + 1e-6
            # Round lightly to keep deterministic, non-integer endpoints
            out.append((round(L, 6), round(R, 6)))
        return out

    # ------------------------ Greedy mex scheduler ----------------------------
    def greedy_mex_order(intervals):
        """
        Adversarial arrival order:
        1) Prime with a short nested chain at a peak-overlap point.
        2) Iteratively pick the interval whose current FirstFit color (mex of blocked)
           is maximal; update blocked colors for neighbors via bit-sets.
        """
        n = len(intervals)

        # Build overlap graph
        neighbors = [[] for _ in range(n)]
        for i in range(n):
            li, ri = intervals[i]
            for j in range(i + 1, n):
                lj, rj = intervals[j]
                if (li < rj) and (lj < ri):
                    neighbors[i].append(j)
                    neighbors[j].append(i)

        # Clique peak point (for priming)
        _, peak_x = clique_size(intervals)
        tower = [i for i, (l, r) in enumerate(intervals) if l < peak_x < r]
        # Sort by length to produce a short nested chain first
        tower.sort(key=lambda idx: (intervals[idx][1] - intervals[idx][0], intervals[idx][0], intervals[idx][1]))
        prime_len = min(5, len(tower))  # quick base palette; modest to avoid stressing omega hot spots

        # Bitset of blocked colors for each node; colors are 1-based -> use bit position "color"
        blocked = [0] * n
        placed = []
        placed_flag = [False] * n
        color_of = [0] * n  # 0 means uncolored/not placed

        def mex(bits):
            # Small colors dominate; typical values under a few dozen
            c = 1
            while bits & (1 << c):
                c += 1
            return c

        # Prime step: place a short chain at the peak
        for idx in tower[:prime_len]:
            # If already placed (unlikely in this loop), skip
            if placed_flag[idx]:
                continue
            c = mex(blocked[idx])
            placed.append(idx)
            placed_flag[idx] = True
            color_of[idx] = c
            # Update neighbors' blocked sets
            for nb in neighbors[idx]:
                if not placed_flag[nb]:
                    blocked[nb] |= (1 << c)

        # Main greedy mex phase
        remaining = n - len(placed)
        # Precompute candidate list for speed
        candidates = [i for i in range(n) if not placed_flag[i]]

        while remaining > 0:
            # Pick interval with maximal mex(blocked)
            best_idx = None
            best_score = -1
            # Tie-breakers: shorter length first, then smaller left endpoint
            best_key = None

            for j in candidates:
                if placed_flag[j]:
                    continue
                score = mex(blocked[j])
                if score > best_score:
                    best_score = score
                    L = intervals[j][0]
                    length = intervals[j][1] - intervals[j][0]
                    best_key = (length, L)
                    best_idx = j
                elif score == best_score:
                    L = intervals[j][0]
                    length = intervals[j][1] - intervals[j][0]
                    key = (length, L)
                    if key < best_key:
                        best_key = key
                        best_idx = j

            # Place the chosen interval
            c = mex(blocked[best_idx])
            placed.append(best_idx)
            placed_flag[best_idx] = True
            color_of[best_idx] = c
            # Update neighbors
            for nb in neighbors[best_idx]:
                if not placed_flag[nb]:
                    blocked[nb] |= (1 << c)

            remaining -= 1

        # Return intervals in placed order
        return [intervals[i] for i in placed]

    # ------------------------ Build, normalize, order -------------------------
    backbone = build_backbone(depth=4)
    normalized = normalize(backbone, target_span=2000.0)
    ordered = greedy_mex_order(normalized)
    return ordered

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()