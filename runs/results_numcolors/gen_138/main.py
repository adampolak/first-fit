# EVOLVE-BLOCK-START

def _rounds_allowed_under_cap(max_intervals, base_size=1, copies=4, connectors=4, target=6):
    """
    Determine how many rounds of the recurrence size_{k+1}=copies*size_k+connectors
    can fit under max_intervals, up to `target` rounds.
    """
    size = base_size
    rounds = 0
    while rounds < target:
        nxt = size * copies + connectors
        if nxt > max_intervals:
            break
        size = nxt
        rounds += 1
    return rounds

def _build_blocks(T, delta, lo, starts, reverse_block_parity):
    """
    Build translated copies of T at offsets given by `starts`.
    Reverse every other block if reverse_block_parity=True.
    """
    blocks = []
    for idx, s in enumerate(starts):
        src = T[::-1] if (reverse_block_parity and (idx % 2 == 1)) else T
        shift = s * delta - lo
        blocks.append([(l + shift, r + shift) for (l, r) in src])
    return blocks

def _interleave_blocks(blocks, interleave, round_idx):
    """
    If interleave=False, flatten blocks in order.
    Otherwise, round-robin across blocks, rotating visitation order each round.
    """
    if not interleave:
        S = []
        for blk in blocks:
            S.extend(blk)
        return S

    maxlen = max((len(b) for b in blocks), default=0)
    order = list(range(len(blocks)))
    # reverse order on odd rounds
    if round_idx % 2 == 1:
        order.reverse()
    # cyclic rotate the order to further break patterns
    k = round_idx % len(order)
    order = order[k:] + order[:k]

    S = []
    for i in range(maxlen):
        for idx in order:
            blk = blocks[idx]
            if i < len(blk):
                S.append(blk[i])
    return S

def _add_connectors(S, starts, delta):
    """
    Add the classic four Kierstead–Trotter connectors (caps) derived
    from the first four values of `starts`.
    """
    s0, s1, s2, s3 = starts[:4]
    connectors = [
        ((s0 - 1) * delta, (s1 - 1) * delta),
        ((s2 + 2) * delta, (s3 + 2) * delta),
        ((s0 + 2) * delta, (s2 - 1) * delta),
        ((s1 + 2) * delta, (s3 - 1) * delta),
    ]
    S.extend(connectors)

def _append_micro_phase(T, delta, max_extra=8):
    """
    Append a small micro-phase of sparse long-range caps near the end.
    """
    if max_extra <= 0:
        return
    lo = min(l for l, _ in T)
    hi = max(r for _, r in T)
    d2 = max(1, delta // 4)
    micro = [
        (lo + 1*d2, lo + 5*d2),
        (hi - 6*d2, hi - 2*d2),
        (lo + 3*d2, lo + 8*d2),
        (hi - 8*d2, hi - 3*d2),
    ]
    for i, iv in enumerate(micro):
        pos = len(T) - (i*2 + 1)
        if pos < 0:
            T.append(iv)
        else:
            T.insert(pos, iv)

def construct_intervals(max_intervals=9800,
                        depth=6,
                        rotate_starts=True,
                        reverse_block_parity=True,
                        interleave_blocks=True,
                        micro_phase=True):
    """
    Crossover Kierstead–Trotter expansion:
      - Adaptive depth under max_intervals.
      - Rotating start patterns to diversify geometry.
      - Optional block-halving reversal, interleaving, in-round gadget, and end micro-phase.
    Returns a list of integer-open intervals [(l,r),...].
    """
    # A small cycle of four-start patterns to rotate each round
    start_patterns = [
        [2, 6, 10, 14],
        [1, 5, 9, 13],
        [3, 7, 11, 15],
        [2, 4, 8, 12],
    ]

    # seed with a single unit interval to keep initial omega=1
    T = [(0.0, 1.0)]

    # figure out how many full rounds we can afford
    rounds = _rounds_allowed_under_cap(max_intervals, base_size=1,
                                       copies=4, connectors=4,
                                       target=depth)

    for r in range(rounds):
        lo = min(l for l, _ in T)
        hi = max(r for _, r in T)
        delta = hi - lo if hi > lo else 1.0

        # pick start pattern this round
        starts = start_patterns[r % len(start_patterns)] if rotate_starts else start_patterns[0]

        # build the four translated blocks
        blocks = _build_blocks(T, delta, lo, starts, reverse_block_parity)

        # interleave or concatenate
        S = _interleave_blocks(blocks, interleave_blocks, r)

        # connectors
        _add_connectors(S, starts, delta)

        # inject a small in-round gadget on the last expansion to boost FF early
        if micro_phase and r == rounds - 1:
            # three fraction-based caps
            for a_frac, b_frac in [(0.08, 0.60), (0.25, 0.75), (0.75, 0.92)]:
                l0 = lo + a_frac * delta
                r0 = lo + b_frac * delta
                if r0 > l0:
                    S.append((l0, r0))

        T = S

    # final micro-phase sprinkling
    if micro_phase and len(T) + 8 <= max_intervals:
        # reuse last delta
        lo = min(l for l, _ in T)
        hi = max(r for _, r in T)
        delta = max(1.0, hi - lo)
        _append_micro_phase(T, delta, max_extra=8)

    # normalize to nonnegative integer endpoints
    if not T:
        return []
    min_l = min(l for l, _ in T)
    if min_l < 0:
        T = [(l - min_l, r - min_l) for l, r in T]

    out = []
    for (l, r) in T:
        li = int(round(l))
        ri = int(round(r))
        if ri <= li:
            ri = li + 1
        out.append((li, ri))

    # truncate to budget
    if len(out) > max_intervals:
        out = out[:max_intervals]
    return out

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()