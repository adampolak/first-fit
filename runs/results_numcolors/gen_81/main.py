# EVOLVE-BLOCK-START

def construct_intervals(max_intervals=9800,
                        depth=6,
                        rotate_starts=True,
                        reverse_block_parity=True,
                        interleave_blocks=True,
                        micro_phase=True):
    """
    Cross-over of rotated Kierstead–Trotter expansions.
    Returns a list of (l, r) integer tuples (open intervals) in presentation order.
    """
    # Rotating start-patterns to vary block offsets
    start_patterns = [
        [2, 6, 10, 14],
        [1, 5,  9, 13],
        [3, 7, 11, 15],
        [2, 4,  8, 12],
    ]
    # Figure 4 gadget to inject in each round
    gadget = [(1, 5), (12, 16), (4, 9), (8, 13)]

    # Seed with a single unit interval
    T = [(0, 1)]

    # Adaptively cap depth so total intervals ≤ max_intervals
    size = len(T)
    allowed = 0
    while allowed < depth:
        # each round adds 4 blocks of size 'size' plus len(gadget) gadgets plus 4 connectors
        next_size = size * 4 + len(gadget) + 4
        if next_size > max_intervals:
            break
        size = next_size
        allowed += 1
    depth = allowed

    for round_idx in range(depth):
        lo = min(l for l, r in T)
        hi = max(r for l, r in T)
        delta = hi - lo or 1

        # choose start-pattern
        if rotate_starts:
            starts = start_patterns[round_idx % len(start_patterns)]
        else:
            starts = start_patterns[0]

        # build four translated blocks
        blocks = []
        for i, s in enumerate(starts):
            src = T[::-1] if (reverse_block_parity and (i % 2 == 1)) else T
            base = s * delta - lo
            blocks.append([(l + base, r + base) for (l, r) in src])

        # interleave or append sequentially
        S = []
        if interleave_blocks:
            maxlen = max(len(b) for b in blocks)
            for i in range(maxlen):
                # rotate block order each round
                for j in range(len(blocks)):
                    blk = blocks[(j + round_idx) % len(blocks)]
                    if i < len(blk):
                        S.append(blk[i])
        else:
            for blk in blocks:
                S.extend(blk)

        # inject Figure 4 gadget scaled by delta
        for a, b in gadget:
            S.append((delta * a, delta * b))

        # connectors (preserve pressure without enlarging ω)
        s0, s1, s2, s3 = starts[:4]
        connectors = [
            ((s0 - 1) * delta, (s1 - 1) * delta),
            ((s2 + 2) * delta, (s3 + 2) * delta),
            ((s0 + 2) * delta, (s2 - 1) * delta),
            ((s1 + 2) * delta, (s3 - 1) * delta),
        ]
        S.extend(connectors)

        T = S

    # micro-phase: sprinkle a few sparse long-range caps
    if micro_phase and len(T) + 4 <= max_intervals:
        lo = min(l for l, r in T)
        hi = max(r for l, r in T)
        delta = hi - lo or 1
        d2 = max(1, delta // 4)
        micro = [
            (lo + 1 * d2, lo + 5 * d2),
            (hi - 6 * d2, hi - 2 * d2),
            (lo + 3 * d2, lo + 8 * d2),
            (hi - 8 * d2, hi - 3 * d2),
        ]
        T.extend(micro)

    # normalize to nonnegative integer intervals
    if T:
        min_l = min(l for l, r in T)
        if min_l < 0:
            T = [(l - min_l, r - min_l) for l, r in T]

    intervals = []
    for l, r in T:
        li = int(round(l))
        ri = int(round(r))
        if ri <= li:
            ri = li + 1
        intervals.append((li, ri))
        if len(intervals) >= max_intervals:
            break

    return intervals

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()