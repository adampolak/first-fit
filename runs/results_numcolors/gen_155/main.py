# EVOLVE-BLOCK-START

def construct_intervals(rounds=6,
                        rotate_starts=True,
                        reverse_block_parity=True,
                        interleave_blocks=True,
                        phase2_iters=2):
    """
    Alternating‐density KT spine with cross4 connectors and a two‐round micro‐phase.
    Same signature and output as the original construct_intervals.
    """
    CAP = 9800

    # Four‐template KT bank
    template_bank = [
        (2, 6, 10, 14),
        (1, 5, 9, 13),
        (3, 7, 11, 15),
        (4, 8, 12, 16),
    ]
    default_starts = template_bank[0]

    # Seed: single unit interval
    T = [(0, 1)]

    # Utility: span and delta of a set
    def _span_delta(S):
        lo = min(l for l, r in S)
        hi = max(r for l, r in S)
        d = hi - lo
        return lo, hi, max(1, d)

    # Build one spine round with multiplier K and optional cross4
    def apply_spine_round(S_prev, starts, do_inter, rev_order, K, add_cross4):
        lo, hi, delta0 = _span_delta(S_prev)
        delta = delta0 * K
        # build translated blocks
        blocks = []
        for s in starts:
            base = s * delta - lo
            block = [(l + base, r + base) for (l, r) in S_prev]
            blocks.append(block)
        # interleave or sequential
        order = list(range(len(blocks)))
        if rev_order:
            order.reverse()
        newS = []
        if do_inter:
            maxlen = max(len(b) for b in blocks)
            for i in range(maxlen):
                for idx in order:
                    blk = blocks[idx]
                    if i < len(blk):
                        newS.append(blk[i])
        else:
            for idx in order:
                newS.extend(blocks[idx])
        # classic four connectors
        s0, s1, s2, s3 = starts
        newS.append(((s0 - 1) * delta, (s1 - 1) * delta))
        newS.append(((s2 + 2) * delta, (s3 + 2) * delta))
        newS.append(((s0 + 2) * delta, (s2 - 1) * delta))
        newS.append(((s1 + 2) * delta, (s3 - 1) * delta))
        # optional cross4
        if add_cross4:
            newS.append(((s0 + 4) * delta, (s3 + 4) * delta))
        return newS

    # Stage 1: spine
    lo0, hi0, _ = _span_delta(T)
    # compute max rounds that stay under CAP
    size = len(T)
    depth = 0
    for _ in range(rounds):
        if 4 * size + 5 > CAP:
            break
        size = 4 * size + 5
        depth += 1

    for ridx in range(depth):
        starts = template_bank[ridx % 4] if rotate_starts else default_starts
        do_inter = interleave_blocks and (ridx % 2 == 0)
        rev_order = reverse_block_parity and (ridx % 2 == 1)
        K = 1 + (ridx % 2)  # alternate density
        add_cross4 = (ridx % 2 == 1) or (ridx == depth - 1)
        T = apply_spine_round(T, starts, do_inter, rev_order, K, add_cross4)
        if len(T) >= CAP:
            T = T[:CAP]
            return T

    # Stage 2: micro‐phase with two delta₂ rounds
    def build_micro_round(S_prev, iter_id, budget):
        lo, hi, D = _span_delta(S_prev)
        delta2 = max(1, D // 3)
        # thin seed
        step = max(1, len(S_prev) // 50)
        U = S_prev[::step][:32]
        if not U:
            return []
        ulo = min(l for l, r in U)
        # choose micro starts patterns
        micro_patterns = [(1, 4), (2, 5)]
        starts = micro_patterns[iter_id % len(micro_patterns)]
        # build micro blocks
        blocks = []
        for idx, s in enumerate(starts):
            base = s * delta2 - ulo
            blk = [(l + base, r + base) for (l, r) in U]
            # small internal parity reversal
            if (iter_id + idx) % 2 == 1:
                blk.reverse()
            blocks.append(blk)
        # interleave parity
        order = list(range(len(blocks)))
        if iter_id % 2 == 1:
            order.reverse()
        M = []
        maxlen = max(len(b) for b in blocks)
        for i in range(maxlen):
            for idx in order:
                if i < len(blocks[idx]):
                    M.append(blocks[idx][i])
        # connectors + cross4 at micro scale
        s0, s1 = starts[0], starts[1]
        M.append((lo + (s0 - 1) * delta2, lo + (s1 - 1) * delta2))
        M.append((lo + (s1 + 2) * delta2, lo + (s0 + 2) * delta2))
        # micro cross4
        M.append((lo + (s0 + 4) * delta2, lo + (s1 + 4) * delta2))
        # enforce budget
        return M[:budget]

    for it in range(min(phase2_iters, 2)):
        room = CAP - len(T)
        if room <= 10:
            break
        micro = build_micro_round(T, it, room)
        if not micro:
            break
        T.extend(micro)

    # final trim
    return T[:CAP]

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()