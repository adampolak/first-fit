# EVOLVE-BLOCK-START

def construct_intervals(seed_count=1):
    """
    Mixed‐arity spine + adaptive micro‐phases.
    Alternates 5‐ and 4‐block KT‐style expansions to diversify FF pressure,
    then runs two delta‐scale micro‐phases tuned to pack new colors without
    blowing up the global clique.
    """

    CAP = 9800

    # Round arities: alternate between richer 5‐block and classic 4‐block expansions
    arities = [5, 4, 5, 4, 5, 4]

    # Seed initialization
    if seed_count <= 1:
        T = [(0, 1)]
    else:
        # small multi‐seed sparsed out to avoid early clique blowup
        cnt = min(seed_count, 5)
        T = [(i * 3, i * 3 + 1) for i in range(cnt)]

    def _span(T):
        lo = min(l for l, r in T)
        hi = max(r for l, r in T)
        d = hi - lo
        return lo, hi, d if d > 0 else 1

    def _apply_round(T, K, do_interleave, reverse_order):
        lo, hi, d = _span(T)
        blocks = []
        # build K translated copies
        for i in range(K):
            # small stagger for odd blocks to break symmetry
            stagger = ((i % 2) * (d // (2 * K)))  
            base = i * d + stagger - lo
            blocks.append([(l + base, r + base) for (l, r) in T])

        S = []
        if do_interleave:
            order = list(range(K))
            if reverse_order:
                order.reverse()
            m = max(len(b) for b in blocks)
            for j in range(m):
                for idx in order:
                    if j < len(blocks[idx]):
                        S.append(blocks[idx][j])
        else:
            if reverse_order:
                blocks.reverse()
            for blk in blocks:
                S.extend(blk)

        # connectors between consecutive blocks to couple colors
        for i in range(K - 1):
            a = lo + i * d + d // 3
            b = lo + (i + 1) * d - d // 3
            if b > a:
                S.append((a, b))

        return S

    # Stage 1: mixed-arity spine
    for ridx, K in enumerate(arities):
        # predictive cap guard
        if 4 * len(T) + K > CAP:
            break
        do_inter = (ridx % 2 == 0)
        rev = (ridx % 3 == 1)
        T = _apply_round(T, K, do_inter, rev)
        if len(T) >= CAP:
            break

    if len(T) >= CAP - 50:
        return T[:CAP]

    # Stage 2: two adaptive micro-phases
    def _micro(T, budget, phase_id):
        lo, hi, d = _span(T)
        n = len(T)
        # choose a thin seed
        seed_sz = min(max(8, n // 200), 50)
        step = max(1, n // seed_sz)
        U = T[::step][:seed_sz]
        if not U or budget < seed_sz:
            return []

        ulo = min(l for l, r in U)
        # build fractional windows (vary count slightly with phase)
        wcount = 4 + (phase_id % 2)
        windows = [(i / (wcount + 1.5), (i + 0.8) / (wcount + 1.5)) for i in range(wcount)]
        blocks = []
        for i, (fa, fb) in enumerate(windows):
            win_lo = lo + int(round(fa * d))
            base = win_lo - ulo
            blk = [(l + base, r + base) for (l, r) in U]
            # occasional reversal to break symmetry
            if (phase_id + i) % 2 == 1:
                blk.reverse()
            blocks.append(blk)

        M = []
        order = list(range(len(blocks)))
        if phase_id % 2 == 1:
            order.reverse()
        m = max(len(b) for b in blocks)
        for j in range(m):
            for idx in order:
                if j < len(blocks[idx]):
                    M.append(blocks[idx][j])

        # window‐scale connectors to force new colors
        for fa, fb in [(0.1, 0.3), (0.4, 0.7), (0.6, 0.9)]:
            a = lo + int(round(fa * d))
            b = lo + int(round(fb * d))
            if b > a:
                M.append((a, b))

        return M[:budget]

    # execute micro-phases
    for phase in range(2):
        if len(T) >= CAP - 10:
            break
        room = CAP - len(T)
        m = _micro(T, room, phase)
        T.extend(m)

    return T[:CAP]

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()