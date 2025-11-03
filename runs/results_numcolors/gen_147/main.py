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

def _append_connectors(starts, delta, base_offset=0, add_cross4=False):
    """Return classic four connectors (caps) and optional long-range cross4."""
    s0, s1, s2, s3 = starts[0], starts[1], starts[2], starts[3]
    connectors = [
        ((s0 - 1) * delta + base_offset, (s1 - 1) * delta + base_offset),
        ((s2 + 2) * delta + base_offset, (s3 + 2) * delta + base_offset),
        ((s0 + 2) * delta + base_offset, (s2 - 1) * delta + base_offset),
        ((s1 + 2) * delta + base_offset, (s3 - 1) * delta + base_offset),
    ]
    if add_cross4:
        connectors.append(((s0 + 4) * delta + base_offset, (s3 + 4) * delta + base_offset))
    return connectors

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

def _make_blocks(T, starts, delta, lo, reverse_parity, K=1):
    """Construct four translated blocks (with optional reversal) and an optional multiplier K."""
    blocks = []
    for b_idx, s in enumerate(starts):
        block_src = T[::-1] if (reverse_parity and (b_idx % 2 == 1)) else T
        base = s * delta * K - lo
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

    if not blocks:
        return []

    maxlen = max((len(b) for b in blocks), default=0)
    order = list(range(len(blocks)))
    if round_idx % 2 == 1:
        order = order[::-1]
    # a subtle rotation to avoid stable overlap patterns
    if len(order) > 0:
        krot = round_idx % len(order)
        order = order[krot:] + order[:krot]

    S = []
    for i in range(maxlen):
        for idx in order:
            blk = blocks[idx]
            if i < len(blk):
                S.append(blk[i])
    return S

def construct_intervals(depth=6, spine_multiplier=1, two_seed=False, add_cross4=True):
    """
    Fractal KT-inspired construction with targeted improvements:
      - spine_multiplier K controls spacing density between translated blocks,
        enabling denser cross-block coupling when K=2 (helps FirstFit).
      - optional two_seed to introduce a small deterministic multi-seed baseline.
      - optional add_cross4 to append a long-range cross4 connector on the final spine round.

    Other heuristics:
      - alternate between interleaved and sequential rounds to both mix active colors
        and create long towers (parity-based interleaving).
      - reverse_block_parity set True to break inner ordering symmetry across blocks.

    Returns: list of (l, r) integer tuples representing open intervals.
    """
    # Global controls (deterministic toggles)
    rotate_starts = True
    reverse_block_parity = True
    # Interleave on even rounds, sequential on odd; this hybrid tends to raise FF
    interleave_blocks = True
    micro_phase = True

    # Start patterns (deterministic cycle, four strong templates)
    start_patterns = [
        [2, 6, 10, 14],
        [1, 5, 9, 13],
        [3, 7, 11, 15],
        [4, 8, 12, 16],
    ]

    # Seed
    if two_seed:
        # Two well-separated unit seeds to encourage early coupling but keep omega low
        T = [(0, 1), (3, 4)]
    else:
        T = [(0, 1)]

    MAX_INTERVALS = 9800  # safety cap

    # Compute a safe depth given the MAX_INTERVALS budget
    depth = max(0, int(depth))
    allowed = 0
    size = len(T)  # current number of intervals notionally
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

        blocks = _make_blocks(T, starts, delta, lo, reverse_block_parity, K=spine_multiplier)

        # Hybrid interleaving: interleave on even rounds, sequential otherwise.
        do_interleave = interleave_blocks and (round_idx % 2 == 0)
        S = _interleave_blocks(blocks, round_idx, interleave=do_interleave)

        # Classic connectors plus an optional long-range cross4 on the final round
        connectors = _append_connectors(starts, delta * spine_multiplier, base_offset=0,
                                        add_cross4=(add_cross4 and round_idx == depth - 1))
        S.extend(connectors)
        T = S

    # Optional micro-phase sprinkled at end (kept small and capacity-bounded)
    if micro_phase and len(T) + 8 <= MAX_INTERVALS:
        lo = min(l for l, _ in T)
        hi = max(r for _, r in T)
        delta = max(1, hi - lo)
        T = _append_micro_phase(T, lo, hi, delta, max_extra=8)

    # EXTRA PHASE: make small secondary pass if we didn't reach target depth
    EXTRA_ROUNDS = 0 if depth >= 6 else min(2, 6 - depth)
    if EXTRA_ROUNDS > 0:
        for extra in range(EXTRA_ROUNDS):
            lo = min(l for l, _ in T)
            hi = max(r for _, r in T)
            delta = hi - lo if hi > lo else 1
            rot = (depth + extra) % len(start_patterns)
            starts = start_patterns[rot]
            blocks = _make_blocks(T, starts, delta, lo, reverse_block_parity, K=spine_multiplier)
            do_interleave = interleave_blocks and ((depth + extra) % 2 == 0)
            S = _interleave_blocks(blocks, depth + extra, interleave=do_interleave)
            connectors = _append_connectors(starts, delta * spine_multiplier, base_offset=0, add_cross4=False)
            S.extend(connectors)
            T = S
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