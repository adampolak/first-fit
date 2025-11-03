# EVOLVE-BLOCK-START

def _rounds_allowed_under_cap(max_intervals, base_size=4, per_round_copies=4, per_round_connectors=4, target_rounds=6):
    """
    Compute how many rounds of the recurrence:
       size_{k+1} = per_round_copies * size_k + per_round_connectors
    fit under max_intervals, up to target_rounds.
    Returns (rounds_used, final_size).
    base_size counts the initial seed block count.
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

def _build_translated_blocks(T, delta, lo, starts):
    """
    Build translated blocks sequentially (no interleaving).
    Each block is a translated copy of T offset by base = s*delta - lo.
    Returns a list of blocks (each block is a list of (l,r)).
    """
    blocks = []
    for s in starts:
        base = s * delta - lo
        blocks.append([(l + base, r + base) for (l, r) in T])
    return blocks

def _append_connectors(S, delta, starts):
    """
    Add the four classic connectors that couple colors across blocks without
    blowing up the clique number. Keeps exactly 4 connectors to fit budgets.
    """
    s0, s1, s2, s3 = (starts[0], starts[1], starts[2], starts[3])
    connectors = [
        ((s0 - 1) * delta, (s1 - 1) * delta),  # left cap
        ((s2 + 2) * delta, (s3 + 2) * delta),  # right cap
        ((s0 + 2) * delta, (s2 - 1) * delta),  # cross 1
        ((s1 + 2) * delta, (s3 - 1) * delta),  # cross 2
    ]
    S.extend(connectors)

def _normalize_to_int_intervals(intervals, cap=None):
    """
    Normalize a list of intervals to non-negative integer endpoints with
    positive length. Optionally truncate to cap elements.
    """
    if not intervals:
        return []
    min_l = min(l for (l, _) in intervals)
    if min_l < 0:
        intervals = [(l - min_l, r - min_l) for (l, r) in intervals]
    out = []
    for (l, r) in intervals:
        li = int(round(l))
        ri = int(round(r))
        if ri <= li:
            ri = li + 1
        out.append((li, ri))
    if cap is not None and len(out) > cap:
        out = out[:cap]
    return out

def _append_micro_phase(T, lo, hi, delta, max_extra=8):
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
    for i, interval in enumerate(micro):
        pos = len(T) - (i * 2 + 1)
        if pos < 0:
            T.append(interval)
        else:
            T.insert(pos, interval)
    return T

def construct_intervals(max_intervals=9800,
                        depth=6,
                        rotate_starts=True,
                        interleave_blocks=True,
                        micro_phase=True):
    """
    Adaptive six-round construction with 4-start blocks, 4 connectors per round,
    and optional interleaving. Seeds with four disjoint unit intervals to boost FF pressure
    while controlling omega through cap budgeting.

    Returns:
      intervals: list[(l, r)] of open intervals with integer endpoints, in arrival order.
    """
    # Deterministic four-pattern cycle (ascending starts)
    start_patterns = [
        (2, 6, 10, 14),  # classic 4-wave
        (1, 5, 9, 13),   # left-shifted
        (3, 7, 11, 15),  # right-shifted
        (2, 4, 8, 12),   # compressed variant
    ]

    # Seed: four disjoint unit intervals
    T = [(0, 1), (2, 3), (4, 5), (6, 7)]
    base_size = len(T)

    # Compute allowed rounds under budget
    rounds_allowed, _ = _rounds_allowed_under_cap(
        max_intervals=max_intervals,
        base_size=base_size,
        per_round_copies=4,
        per_round_connectors=4,
        target_rounds=int(depth)
    )
    if rounds_allowed <= 0:
        return [(0, 1)]

    # Main expansion
    for round_idx in range(rounds_allowed):
        lo = min(l for l, _ in T)
        hi = max(r for _, r in T)
        delta = hi - lo
        if delta <= 0:
            delta = 1

        starts = start_patterns[round_idx % len(start_patterns)] if rotate_starts else start_patterns[0]

        # Build four translated blocks in sequential order
        blocks = _build_translated_blocks(T, delta, lo, starts)

        # Interleave blocks if requested
        S = []
        if interleave_blocks:
            maxlen = max((len(b) for b in blocks), default=0)
            order = list(range(len(blocks)))
            if round_idx % 2 == 1:
                order = order[::-1]
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

        # Deterministic connectors (caps) derived from the first four starts
        _append_connectors(S, delta, starts)

        T = S

    # Optional micro-phase sprinkled at end
    if micro_phase and len(T) + 8 <= max_intervals:
        lo = min(l for l, _ in T)
        hi = max(r for _, r in T)
        delta = max(1, hi - lo)
        T = _append_micro_phase(T, lo, hi, delta, max_extra=8)

    # Normalize and cap
    intervals = _normalize_to_int_intervals(T, cap=max_intervals)
    return intervals

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()