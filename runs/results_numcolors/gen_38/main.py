# EVOLVE-BLOCK-START

def construct_intervals(max_intervals=9800,
                        requested_rounds=6,
                        seed_units=4,
                        rotate_starts=True,
                        reverse_block_parity=True,
                        interleave_blocks=True,
                        micro_phase=True):
    """
    Rotated tower-and-cap style deterministic construction.

    Returns a list of (l, r) integer tuples (open intervals) in the order
    they should be presented to FirstFit.

    Parameters are chosen so that the total number of intervals stays < max_intervals.
    """
    # Deterministic bank of start-pattern templates (4-block offsets)
    templates = [
        [2, 6, 10, 14],   # classic
        [1, 5, 9, 13],    # left shifted
        [3, 7, 11, 15],   # right shifted
        [2, 4, 8, 12],    # compressed variant
    ]

    # Seed: a small richer set of disjoint unit intervals (spaced out)
    spacing = 4
    T = [(i * spacing, i * spacing + 1) for i in range(max(1, int(seed_units)))]

    # Safety: estimate growth and reduce rounds if needed.
    # Use conservative connector count per round to ensure we don't exceed max_intervals.
    est_size = len(T)
    allowed_rounds = 0
    est_rounds = max(0, int(requested_rounds))
    # conservative connector count (we insert several connectors each round)
    connectors_per_round = 12
    while allowed_rounds < est_rounds:
        next_size = est_size * 4 + connectors_per_round
        if next_size > max_intervals:
            break
        est_size = next_size
        allowed_rounds += 1

    rounds = allowed_rounds

    # Main rotated expansion
    for round_idx in range(rounds):
        lo = min(l for l, r in T)
        hi = max(r for l, r in T)
        delta = hi - lo
        if delta <= 0:
            delta = 1

        starts = templates[round_idx % len(templates)] if rotate_starts else templates[0]

        # Build translated copies (possibly reversing every other block)
        blocks = []
        for b_idx, s in enumerate(starts):
            src = T[::-1] if (reverse_block_parity and (b_idx % 2 == 1)) else T
            base = s * delta - lo
            block = [(l + base, r + base) for (l, r) in src]
            blocks.append(block)

        # Interleave blocks (round-robin) to maximize mixing of colors
        S = []
        if interleave_blocks:
            maxlen = max(len(b) for b in blocks)
            order = list(range(len(blocks)))
            # vary the visitation order deterministically
            if round_idx % 2 == 1:
                order = order[::-1]
            # rotate the order each round for additional disruption
            krot = round_idx % len(order)
            order = order[krot:] + order[:krot]
            for i in range(maxlen):
                for idx in order:
                    blk = blocks[idx]
                    if i < len(blk):
                        S.append(blk[i])
        else:
            for blk in blocks:
                S.extend(blk)

        # Build a diverse set of connectors (left/right/cross/adjacent/bridges).
        # Use start positions to compute geometrically-placed connectors.
        s0, s1, s2, s3 = starts
        raw_connectors = [
            ((s0 - 1) * delta, (s1 - 1) * delta),  # left cap
            ((s2 + 2) * delta, (s3 + 2) * delta),  # right cap
            ((s0 + 2) * delta, (s2 - 1) * delta),  # cross 1
            ((s1 + 2) * delta, (s3 - 1) * delta),  # cross 2
        ]
        # Adjacent bridges between neighboring translated copies
        for i in range(3):
            a = (starts[i] + 0.5) * delta
            b = (starts[i + 1] - 0.5) * delta
            raw_connectors.append((a, b))
        # Add a couple extra long-range couplers to propagate colors further
        raw_connectors.append(((starts[0] - 0.5) * delta, (starts[3] + 0.5) * delta))
        raw_connectors.append(((starts[1] - 0.25) * delta, (starts[2] + 0.25) * delta))

        # Normalize connectors to valid intervals with positive length
        connectors = []
        for a, b in raw_connectors:
            l_conn = float(min(a, b))
            r_conn = float(max(a, b))
            # make sure connectors have at least a small positive span
            if r_conn - l_conn < 1.0:
                r_conn = l_conn + 1.0
            connectors.append((l_conn, r_conn))

        # Insert connectors distributed into S (not just appended) to increase mixing.
        if connectors:
            step = max(1, len(S) // (len(connectors) + 1))
            insert_pos = step
            for conn in connectors:
                if insert_pos >= len(S):
                    S.append(conn)
                else:
                    S.insert(insert_pos, conn)
                insert_pos += step + 1

        T = S

    # Micro-phase: add a small-scale second wave of gadgets and a set of sparse, long caps.
    if micro_phase:
        # compute current span
        if len(T) >= 1:
            lo = min(l for l, r in T)
            hi = max(r for l, r in T)
            delta = max(1, hi - lo)
        else:
            lo, hi, delta = 0.0, 1.0, 1.0

        # Choose a small base subset of T to clone on micro-scale (polish few hotspots)
        subset_size = max(4, min(128, len(T) // 32))
        base_subset = T[:subset_size]

        # micro scale factor
        d2 = max(1.0, delta / 5.0)

        # deterministic micro templates (small shifts)
        micro_starts = [1, 3, 5, 7]
        micro_inserts = []
        for m_idx, ms in enumerate(micro_starts):
            shift = lo + ms * d2
            for (l, r) in base_subset:
                # compress the width slightly to create smaller micro-intervals
                mid = (l + r) / 2.0
                half = max(0.3, (r - l) * 0.4)
                nl = mid - half + shift * 0.0001  # tiny deterministic perturbation
                nr = mid + half + shift * 0.0002
                micro_inserts.append((nl, nr))
            # a short connector tying this micro-block into the main structure
            micro_inserts.append((shift - d2 * 0.5, shift + d2 * 1.5))

        # Insert micro gadgets evenly just before the end to maximize pressure late.
        insert_base = max(0, len(T) - len(micro_inserts) - 5)
        for i, it in enumerate(micro_inserts):
            pos = insert_base + i
            if pos >= len(T):
                T.append(it)
            else:
                T.insert(pos, it)

        # Add a few long sparse caps (these intersect many earlier translated copies,
        # but are placed so as not to blow up clique size).
        caps = []
        span = max(1.0, hi - lo)
        # create caps of different lengths covering staggered subranges
        for k in range(12):
            a = lo - span * 0.1 + (k * span * 0.07)
            b = a + span * (0.18 + 0.02 * (k % 3))
            # ensure proper ordering and positive width
            if b - a < 1.0:
                b = a + 1.0
            caps.append((a, b))
        # Append caps near the end
        for cap in caps:
            T.append(cap)

    # Final normalization: shift so min endpoint >= 0 and convert to integers
    if not T:
        return []
    min_l = min(l for l, r in T)
    if min_l < 0:
        shift_amt = -min_l
        T = [(l + shift_amt, r + shift_amt) for (l, r) in T]

    # Convert to integer endpoints while ensuring positive lengths.
    intervals = []
    for (l, r) in T:
        # round endpoints to nearest integer but maintain r > l
        li = int(round(l))
        ri = int(round(r))
        if ri <= li:
            ri = li + 1
        intervals.append((li, ri))

    # Ensure we did not exceed the budget; if we did, truncate some trailing caps/gadgets.
    if len(intervals) > max_intervals:
        intervals = intervals[:max_intervals]

    return intervals


def run_experiment(**kwargs):
    """Main called by evaluator"""
    # allow external override of rounds/seed via kwargs, but keep defaults tuned
    return construct_intervals(**kwargs)

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()