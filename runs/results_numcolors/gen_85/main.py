# EVOLVE-BLOCK-START

def _rounds_allowed_under_cap(max_intervals, base_size=1, per_round_copies=4, per_round_connectors=4, target_rounds=6):
    """
    Compute how many rounds of the recurrence:
       size_{k+1} = per_round_copies * size_k + per_round_connectors
    fit under max_intervals, up to target_rounds.
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
    This ordering is known to increase FF pressure.
    """
    S = []
    for s in starts:
        shift = s * delta - lo
        S.extend((l + shift, r + shift) for (l, r) in T)
    return S

def _append_connectors(S, delta, starts):
    """
    Add the four classic connectors that couple colors across blocks without
    blowing up the clique number. Keep exactly 4 per round to fit the 6-round budget.
    """
    s0, s1, s2, s3 = starts
    connectors = [
        ((s0 - 1) * delta, (s1 - 1) * delta),  # left cap
        ((s2 + 2) * delta, (s3 + 2) * delta),  # right cap
        ((s0 + 2) * delta, (s2 - 1) * delta),  # cross 1
        ((s1 + 2) * delta, (s3 - 1) * delta),  # cross 2
    ]
    S.extend(connectors)

def _normalize_non_negative_ints(intervals, cap=None):
    """
    Shift all intervals to be non-negative and round to integers,
    ensuring open intervals have positive length.
    Optionally truncate to 'cap' total intervals.
    """
    if not intervals:
        return []
    min_l = min(l for l, _ in intervals)
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

def construct_intervals(max_intervals=9800,
                        depth=6,
                        rotate_starts=True,
                        reverse_block_parity=True,   # kept for interface; not used in the classic recipe
                        interleave_blocks=False,     # default to sequential, which empirically boosts FF
                        micro_phase=True):
    """
    Classic four-start, six-round recursive construction with 4 connectors per round.
    This reliably yields a high FF/omega ratio under the 10k-interval budget.

    Parameters:
      max_intervals (int): overall cap on number of intervals.
      depth (int): requested rounds (we will try to reach at least 6 if budget allows).
      rotate_starts (bool): cycle through deterministic four-pattern start templates.
      reverse_block_parity (bool): retained for compatibility (unused).
      interleave_blocks (bool): retained for compatibility (unused; sequential is used).
      micro_phase (bool): optional post-phase small caps (skipped if 6 rounds are used).

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

    # We prefer six rounds if budget allows, regardless of 'depth' (while respecting the cap).
    target_rounds = max(6, int(depth))

    # Compute allowed rounds under the cap for the 4x + 4 recurrence
    rounds_allowed, _ = _rounds_allowed_under_cap(
        max_intervals=max_intervals,
        base_size=1,
        per_round_copies=4,
        per_round_connectors=4,
        target_rounds=target_rounds
    )
    # If for some reason we can't fit 1 round (extremely small cap), bail out early
    if rounds_allowed <= 0:
        return [(0, 1)]

    # Seed: a single unit interval; this keeps omega small and enables many rounds within budget
    T = [(0, 1)]

    # Main expansion
    for round_idx in range(rounds_allowed):
        lo = min(l for l, r in T)
        hi = max(r for l, r in T)
        delta = hi - lo
        if delta <= 0:
            delta = 1

        starts = start_patterns[round_idx % len(start_patterns)] if rotate_starts else start_patterns[0]

        # Build four translated blocks in sequential order (no interleaving)
        S = _build_translated_blocks(T, delta, lo, starts)

        # Append the four deterministic connectors
        _append_connectors(S, delta, starts)

        T = S

    # Optional micro-phase: only if we did not reach 6 rounds
    # (The 6-round construction already achieves a strong ratio; extra caps may disrupt it.)
    if micro_phase and rounds_allowed < 6:
        lo = min(l for l, r in T)
        hi = max(r for l, r in T)
        delta = max(1, hi - lo)
        d2 = max(1, delta // 4)
        micro = [
            (lo + 1 * d2, lo + 5 * d2),
            (hi - 6 * d2, hi - 2 * d2),
            (lo + 3 * d2, lo + 8 * d2),
            (hi - 8 * d2, hi - 3 * d2),
        ]
        # Insert these near the end to touch many active colors
        for i, it in enumerate(micro):
            pos = len(T) - (2 * i + 1)
            if 0 <= pos < len(T):
                T.insert(pos, it)
            else:
                T.append(it)

    # Normalize and truncate if needed
    intervals = _normalize_non_negative_ints(T, cap=max_intervals)
    return intervals

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()