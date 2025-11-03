# EVOLVE-BLOCK-START

def construct_intervals(depth=5,
                        rotate_starts=True,
                        reverse_block_parity=True,
                        interleave_blocks=True,
                        preload_blockers=2,
                        blocker_length_frac=0.03,
                        blocker_gap_frac=0.02):
    """
    Enhanced recursive four-block expansion inspired by Kierstead-Trotter (Figure 4).
    Improvements over the prior version:
      - rotation of start patterns across rounds to disrupt repeating overlaps,
      - optional reversal of block order for odd blocks to mix FirstFit color assignment,
      - small preload blocker intervals inside each translated block to occupy low colors early,
      - short bridges between adjacent blocks to couple colors,
      - default depth increased to 5 (still under the usual 10k cap).
    Parameters are conservative; they can be tuned to push FirstFit further.
    Returns a list of (l,r) intervals (integers) in presentation order.
    """
    # Deterministic rotation cycle for the four translated copies per round
    start_patterns = [
        [2, 6, 10, 14],  # classic
        [1, 5, 9, 13],   # left-shifted
        [3, 7, 11, 15],  # right-shifted
        [2, 4, 8, 12],   # compressed variant
    ]

    # Base seed
    T = [(0.0, 1.0)]

    # Guard depth
    depth = max(0, int(depth))

    for round_idx in range(depth):
        lo = min(l for l, r in T)
        hi = max(r for l, r in T)
        span = hi - lo
        if span <= 0:
            span = 1.0

        starts = start_patterns[round_idx % len(start_patterns)] if rotate_starts else start_patterns[0]

        blocks = []
        # Build each translated block, optionally reversing the clone order
        for b_idx, s in enumerate(starts):
            base = s * span - lo
            block = []
            # insert small preload blockers at the start of the block area to occupy small colors
            for b in range(preload_blockers):
                start = base + (0.06 + b * (blocker_length_frac + blocker_gap_frac)) * span
                end = start + blocker_length_frac * span
                block.append((start, end))
            # append the cloned copy of T (possibly reversed)
            src = list(reversed(T)) if (reverse_block_parity and (b_idx % 2 == 1)) else T
            for (l, r) in src:
                block.append((l + base, r + base))
            blocks.append(block)

        # Optionally interleave blocks round-robin to mix colors more aggressively
        S = []
        if interleave_blocks:
            maxlen = max(len(b) for b in blocks)
            order = list(range(len(blocks)))
            # rotate the order deterministically each round to vary interactions
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

        # Add short bridges between adjacent blocks to couple colors lightly
        for i in range(len(starts) - 1):
            s = starts[i]
            next_s = starts[i + 1]
            S.append(((s + 0.55) * span, (next_s + 0.35) * span))

        # Add the deterministic connector intervals (caps / crosses) derived from the classic pattern
        s0, s1, s2, s3 = starts[0], starts[1], starts[2], starts[3]
        connectors = [
            ((s0 - 1) * span, (s1 - 1) * span),  # left cap
            ((s2 + 2) * span, (s3 + 2) * span),  # right cap
            ((s0 + 2) * span, (s2 - 1) * span),  # cross 1
            ((s1 + 2) * span, (s3 - 1) * span),  # cross 2
        ]
        for a, b in connectors:
            S.append((a, b))

        T = S

    # Convert to integer coordinates with scaling to preserve fractional structure
    if not T:
        return []
    # scale factor chosen to keep relative differences visible when rounding
    scale_factor = 1000
    min_l = min(l for l, r in T)
    # shift to non-negative
    shifted = [( (l - min_l) * scale_factor, (r - min_l) * scale_factor ) for (l, r) in T]

    intervals = []
    for (l, r) in shifted:
        li = int(round(l))
        ri = int(round(r))
        if ri <= li:
            ri = li + 1
        intervals.append((li, ri))

    # keep the sequence size reasonable (safety cap)
    if len(intervals) > 10000:
        intervals = intervals[:10000]

    return intervals

def run_experiment(**kwargs):
    return construct_intervals()

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()