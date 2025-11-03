# EVOLVE-BLOCK-START

def _rounds_allowed_under_cap(max_intervals, base_size=1, per_round_copies=4, per_round_connectors=4, target_rounds=6):
    """
    Compute how many rounds of the recurrence:
       size_{k+1} = per_round_copies * size_k + per_round_connectors
    fit under max_intervals, up to target_rounds.
    Returns (rounds_used, final_size).
    """
    size = base_size
    rounds = 0
    while rounds < target_rounds:
        next_size = per_round_copies * size + per_round_connectors
        if next_size > max_intervals:
            break
        size = next_size
        rounds += 1
    return rounds, size

def _build_rotated_blocks(T, delta, lo, starts, reverse_block_parity):
    """
    Build 4 translated blocks from T in a deterministic rotated pattern.
    Each block is shifted by start*delta - lo, and every other block can be reversed
    to disrupt color reuse and increase FF pressure.
    Returns a list of 4 blocks, each a list of (l, r) pairs.
    """
    blocks = []
    for b_idx, s in enumerate(starts):
        block_src = T[::-1] if (reverse_block_parity and (b_idx % 2 == 1)) else T
        shift = s * delta - lo
        block = [(l + shift, r + shift) for (l, r) in block_src]
        blocks.append(block)
    return blocks

def _interleave_blocks(blocks, round_idx, interleave=True):
    """
    Interleave blocks in a round-robin fashion to maximize cross-block overlap.
    If interleave is False, concatenate blocks in order.
    The visitation order is alternated every round to reduce predictable reuse.
    """
    if not interleave:
        S = []
        for blk in blocks:
            S.extend(blk)
        return S

    maxlen = max((len(b) for b in blocks), default=0)
    order = list(range(len(blocks)))
    if round_idx % 2 == 1:
        order = order[::-1]
    S = []
    for i in range(maxlen):
        for idx in order:
            blk = blocks[idx]
            if i < len(blk):
                S.append(blk[i])
    return S

def _add_connectors(S, starts, delta):
    """
    Add the four deterministic connectors (caps) that couple blocks without
    creating a large clique. The formulas mirror classic Figure-4 style patterns.
    """
    s0, s1, s2, s3 = starts
    connectors = [
        ((s0 - 1) * delta, (s1 - 1) * delta),  # left cap
        ((s2 + 2) * delta, (s3 + 2) * delta),  # right cap
        ((s0 + 2) * delta, (s2 - 1) * delta),  # cross 1
        ((s1 + 2) * delta, (s3 - 1) * delta),  # cross 2
    ]
    S.extend(connectors)

def _normalize_to_non_negative(intervals):
    """
    Shift all intervals so that the minimum left endpoint is non-negative.
    Convert to integer endpoints to ensure clean presentation.
    """
    if not intervals:
        return intervals
    min_l = min(l for (l, _) in intervals)
    if min_l < 0:
        intervals = [(l - min_l, r - min_l) for (l, r) in intervals]
    return [(int(l), int(r)) for (l, r) in intervals]

def construct_intervals(max_intervals=9800,
                      depth=6,
                      rotate_starts=True,
                      reverse_block_parity=True,
                      interleave_blocks=True,
                      phase2_micro=False):
    """
    Rotated four-block Kierstead–Trotter style expansion with deterministic connectors.
    - max_intervals: total allowed intervals (cap).
    - depth: nominal maximum rounds; actual rounds may be smaller due to cap.
    - rotate_starts: rotate the four-start pattern across rounds to vary geometry.
    - reverse_block_parity: reverse every other block to mix color usage.
    - interleave_blocks: interleave blocks round-robin to maximize FF pressure.
    - phase2_micro: optionally sprinkle a tiny micro-phase of long sparse caps.
    Returns:
      intervals: list[(l, r)] of open intervals with integer endpoints, in arrival order.
    """
    # Deterministic four-pattern cycle (ascending starts)
    start_patterns = [
        [2, 6, 10, 14],  # classic
        [1, 5, 9, 13],   # left-shifted
        [3, 7, 11, 15],  # right-shifted
        [2, 4, 8, 12],   # compressed variant
    ]

    # Compute usable rounds under cap
    rounds_allowed, _ = _rounds_allowed_under_cap(
        max_intervals=max_intervals,
        base_size=1,
        per_round_copies=4,
        per_round_connectors=4,
        target_rounds=depth
    )
    if rounds_allowed <= 0:
        return [(0, 1)]

    # Seed: a single unit interval
    T = [(0, 1)]

    for round_idx in range(rounds_allowed):
        lo = min(l for l, _ in T)
        hi = max(r for _, r in T)
        delta = hi - lo
        if delta <= 0:
            delta = 1

        starts = start_patterns[round_idx % len(start_patterns)] if rotate_starts else start_patterns[0]

        # Build four translated blocks in sequential order
        blocks = _build_rotated_blocks(T, delta, lo, starts, reverse_block_parity)

        # Interleave or not
        S = _interleave_blocks(blocks, round_idx, interleave=interleave_blocks)

        # Append connectors
        _add_connectors(S, starts, delta)

        T = S

    # Optional tiny micro-phase
    if phase2_micro and len(T) + 8 <= max_intervals:
        lo = min(l for l, _ in T)
        hi = max(r for _, r in T)
        delta = max(1, hi - lo)
        d2 = max(1, delta // 4)
        micro = [
            (lo + 1 * d2, lo + 5 * d2),
            (hi - 6 * d2, hi - 2 * d2),
            (lo + 3 * d2, lo + 8 * d2),
            (hi - 8 * d2, hi - 3 * d2),
        ]
        # sprinkle near the end, alternating to touch many active colors
        for i, interval in enumerate(micro):
            insert_pos = len(T) - (i * 2 + 1)
            if insert_pos < 0:
                T.append(interval)
            else:
                T.insert(insert_pos, interval)

    intervals = _normalize_to_non_negative(T)
    if len(intervals) > max_intervals:
        intervals = intervals[:max_intervals]
    return intervals

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()