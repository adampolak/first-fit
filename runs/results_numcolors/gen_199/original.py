# EVOLVE-BLOCK-START

def _normalize_to_int_intervals(intervals):
    """Convert to non-negative integer endpoints with minimal adjustments."""
    if not intervals:
        return []
    min_l = min(l for l, _ in intervals)
    if min_l < 0:
        intervals = [(l - min_l, r - min_l) for l, r in intervals]
    # Ensure strictly positive span per interval
    normalized = []
    for (l, r) in intervals:
        li = int(l)
        ri = int(r)
        if ri <= li:
            ri = li + 1
        normalized.append((li, ri))
    return normalized

def _append_micro_phase(T, lo, hi, delta, max_extra):
    """Append a small micro-phase at the end to sprinkle long caps, bounded by max_extra."""
    if max_extra <= 0:
        return T
    d2 = max(1, delta // 4)
    micro = [
        (lo + 1 * d2, lo + 5 * d2),
        (hi - 6 * d2, hi - 2 * d2),
        (lo + 3 * d2, lo + 8 * d2),
        (hi - 8 * d2, hi - 3 * d2),
    ]
    # insert near end to touch many active colors
    for i, interval in enumerate(micro):
        pos = len(T) - (i * 2 + 1)
        if pos < 0:
            T.append(interval)
        else:
            T.insert(pos, interval)
    return T

def _make_blocks(T, starts, delta, lo, reverse_parity):
    """Construct four translated blocks (with optional reversal)"""
    blocks = []
    for b_idx, s in enumerate(starts):
        block_src = T[::-1] if (reverse_parity and (b_idx % 2 == 1)) else T
        base = s * delta - lo
        block = [(l + base, r + base) for (l, r) in block_src]
        blocks.append(block)
    return blocks

def _interleave_blocks(blocks, round_idx, interleave=True):
    """Return an interleaved sequence of intervals from blocks; round-dependent order."""
    if not interleave:
        # flatten in block order
        S = []
        for blk in blocks:
            S.extend(blk)
        return S

    maxlen = max((len(b) for b in blocks), default=0)
    order = list(range(len(blocks)))
    if round_idx % 2 == 1:
        order = order[::-1]
    # a subtle rotation to avoid stable overlap patterns
    krot = round_idx % len(order)
    order = order[krot:] + order[:krot]

    S = []
    for i in range(maxlen):
        for idx in order:
            blk = blocks[idx]
            if i < len(blk):
                S.append(blk[i])
    return S

def _thin_seed(current_T, max_seed):
    """Take a thin, evenly spaced sample of current_T of size <= max_seed."""
    n = len(current_T)
    if n == 0 or max_seed <= 0:
        return []
    step = max(1, n // max_seed)
    U = current_T[::step][:max_seed]
    return U

def _build_micro_delta_round(current_T, budget, iter_id=0, alt=False):
    """
    Build a thin, windowed micro-round translated across the global span.
    Uses two distinct window families across stages; adds deterministic connectors.
    Capacity-guarded by 'budget'.
    """
    if not current_T or budget <= 0:
        return []

    glo = min(l for l, _ in current_T)
    ghi = max(r for _, r in current_T)
    G = max(1, ghi - glo)

    # Thin seed limits omega while keeping strong cross-scale coupling
    seed_sz = max(8, min(48, len(current_T) // 260))
    U = _thin_seed(current_T, seed_sz)
    if not U:
        return []

    ulo = min(l for l, _ in U)

    # Two window families; primary has a tiny shift, alternate is fixed
    if not alt:
        shift = (iter_id % 3) * 0.02
        window_fracs = [
            (0.12 + shift, 0.22 + shift),
            (0.35 + shift, 0.45 + shift),
            (0.58 + shift, 0.68 + shift),
            (0.80 + shift, 0.90 + shift),
        ]
        window_fracs = [(max(0.04, a), min(0.96, b)) for (a, b) in window_fracs]
    else:
        window_fracs = [(0.05, 0.15), (0.28, 0.38), (0.60, 0.70), (0.82, 0.92)]

    # Build translated micro-blocks and parity-dependent internal reversals
    blocks = []
    for (fa, fb) in window_fracs:
        win_lo = glo + int(round(fa * G))
        base = win_lo - ulo
        block = [(l + base, r + base) for (l, r) in U]
        tag = iter_id + (1 if alt else 0)
        if ((int(round(fa * 100)) // 5) + tag) % 2 == 1:
            block = list(reversed(block))
        blocks.append(block)

    # Interleave blocks; reverse order on odd tags
    micro = []
    maxlen = max(len(b) for b in blocks)
    order = list(range(len(blocks)))
    tag = iter_id + (1 if alt else 0)
    if tag % 2 == 1:
        order.reverse()
    for i in range(maxlen):
        for idx in order:
            blk = blocks[idx]
            if i < len(blk):
                micro.append(blk[i])

    # Deterministic connectors across windows (fractional-span analog of KT caps)
    micro_connectors = [
        (glo + int(round(0.08 * G)), glo + int(round(0.30 * G))),  # left cap
        (glo + int(round(0.60 * G)), glo + int(round(0.92 * G))),  # right cap
        (glo + int(round(0.26 * G)), glo + int(round(0.56 * G))),  # cross 1
        (glo + int(round(0.44 * G)), glo + int(round(0.78 * G))),  # cross 2
    ]
    # A longer-range cross only in the alternate micro-phase
    if alt:
        micro_connectors.append((glo + int(round(0.18 * G)), glo + int(round(0.84 * G))))

    for a, b in micro_connectors:
        if b > a:
            micro.append((a, b))

    # Trim to budget
    if len(micro) > budget:
        micro = micro[:budget]
    return micro

def construct_intervals(depth=6):
    """
    Fractal KT-inspired construction with modular toggles.
    Tuned to reach six rounds under the cap and use sequential block order,
    which empirically increases FirstFit pressure without inflating omega.
    Returns: list of (l, r) integer tuples representing open intervals.
    """
    # Global controls (deterministic toggles)
    rotate_starts   = True
    reverse_block_parity = False   # keep inner order stable within blocks
    interleave_blocks = False      # sequential blocks boost FF vs omega
    micro_phase = True

    # Start patterns (deterministic cycle)
    start_patterns = [
        [2, 6, 10, 14],  # classic four-block layout
        [1, 5, 9, 13],   # left-shifted
        [3, 7, 11, 15],  # right-shifted
        [2, 4, 8, 12],   # compressed variant
    ]

    # Seed: a single unit interval; small to allow deeper rounds under cap
    T = [(0, 1)]

    MAX_INTERVALS = 9800  # safety cap

    # Compute a safe depth given the MAX_INTERVALS budget
    depth = max(0, int(depth))
    allowed = 0
    size = 1  # current number of intervals notionally
    while allowed < depth:
        next_size = size * 4 + 4  # 4 copies + 4 connectors (rough bound)
        if next_size > MAX_INTERVALS:
            break
        size = next_size
        allowed += 1
    depth = allowed

    for round_idx in range(depth):
        lo = min(l for l, _ in T)
        hi = max(r for _, r in T)
        delta = hi - lo
        if delta <= 0:
            delta = 1

        starts = start_patterns[round_idx % len(start_patterns)] if rotate_starts else start_patterns[0]

        blocks = _make_blocks(T, starts, delta, lo, reverse_block_parity)

        S = _interleave_blocks(blocks, round_idx, interleave=interleave_blocks)

        # Deterministic four-connectors (caps) derived from the first four starts
        s0, s1, s2, s3 = starts[0], starts[1], starts[2], starts[3]
        connectors = [
            ((s0 - 1) * delta, (s1 - 1) * delta),  # left cap
            ((s2 + 2) * delta, (s3 + 2) * delta),  # right cap
            ((s0 + 2) * delta, (s2 - 1) * delta),  # cross 1
            ((s1 + 2) * delta, (s3 - 1) * delta),  # cross 2
        ]
        S.extend(connectors)

        T = S

    # Optional micro-phase sprinkled at end (tiny long caps near tail)
    if micro_phase and len(T) + 8 <= MAX_INTERVALS:
        lo = min(l for l, _ in T)
        hi = max(r for _, r in T)
        delta = max(1, hi - lo)
        T = _append_micro_phase(T, lo, hi, delta, max_extra=8)

    # Deterministic two-stage delta2 micro-phase to push FF while guarding omega
    room = MAX_INTERVALS - len(T)
    if room > 12:
        budgetA = min(room // 2, 200)
        microA = _build_micro_delta_round(T, budget=budgetA, iter_id=0, alt=False)
        if microA:
            avail = MAX_INTERVALS - len(T)
            if len(microA) > avail:
                microA = microA[:avail]
            T.extend(microA)

    room = MAX_INTERVALS - len(T)
    if room > 12:
        budgetB = min(room, 200)
        microB = _build_micro_delta_round(T, budget=budgetB, iter_id=1, alt=True)
        if microB:
            avail = MAX_INTERVALS - len(T)
            if len(microB) > avail:
                microB = microB[:avail]
            T.extend(microB)

    # EXTRA PHASE: disabled when depth >= 6 (default); fallback if small depths are requested
    EXTRA_ROUNDS = 0 if depth >= 6 else min(2, 6 - depth)
    if EXTRA_ROUNDS > 0:
        for extra in range(EXTRA_ROUNDS):
            lo = min(l for l, _ in T)
            hi = max(r for _, r in T)
            delta = hi - lo
            if delta <= 0:
                delta = 1
            # Rotate starts deterministically for the extra rounds
            rot = (depth + extra) % len(start_patterns)
            starts = start_patterns[rot]

            blocks = _make_blocks(T, starts, delta, lo, reverse_block_parity)
            S = _interleave_blocks(blocks, depth + extra, interleave=interleave_blocks)

            # Deterministic connectors (retain 4 to keep the budget tight)
            s0, s1, s2, s3 = starts[0], starts[1], starts[2], starts[3]
            connectors = [
                ((s0 - 1) * delta, (s1 - 1) * delta),  # left cap
                ((s2 + 2) * delta, (s3 + 2) * delta),  # right cap
                ((s0 + 2) * delta, (s2 - 1) * delta),  # cross 1
                ((s1 + 2) * delta, (s3 - 1) * delta),  # cross 2
            ]
            S.extend(connectors)
            T = S

            # Enforce the interval cap defensively
            if len(T) > MAX_INTERVALS:
                T = T[:MAX_INTERVALS]
                break

    # Final normalization
    intervals = _normalize_to_int_intervals(T)
    return intervals


def run_experiment(**kwargs):
    """Main called by evaluator"""
    return construct_intervals(**kwargs)

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()