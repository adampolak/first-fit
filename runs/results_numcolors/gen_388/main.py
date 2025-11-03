# EVOLVE-BLOCK-START

def construct_intervals(rounds=6,
                        rotate_starts=True,
                        reverse_block_parity=True,
                        interleave_blocks=True,
                        phase2_iters=1,
                        cross4_enabled=True):
    """
    Hybrid KT‐spine plus fractional micro‐phases.
    Returns a list of open intervals (l, r) in FF arrival order.
    """
    CAP = 9800

    # Core templates
    default_starts = (2, 6, 10, 14)
    template_bank = [
        default_starts,
        (1, 5, 9, 13),
        (3, 7, 11, 15),
        (4, 8, 12, 16),
    ]

    # Seed
    T = [(0, 1)]

    # Helpers
    def span_delta(seq):
        lo = min(l for l, r in seq)
        hi = max(r for l, r in seq)
        d = hi - lo
        return lo, hi, d if d > 0 else 1

    def append_connectors(S, starts, delta, add_cross4=False, lo_off=0):
        s0, s1, s2, s3 = starts
        S.append((lo_off + (s0 - 1) * delta, lo_off + (s1 - 1) * delta))
        S.append((lo_off + (s2 + 2) * delta, lo_off + (s3 + 2) * delta))
        S.append((lo_off + (s0 + 2) * delta, lo_off + (s2 - 1) * delta))
        S.append((lo_off + (s1 + 2) * delta, lo_off + (s3 - 1) * delta))
        if add_cross4:
            # sparse long‐range connector
            S.append((lo_off + (s0 + 4) * delta, lo_off + (s3 + 4) * delta))

    def apply_round(cur, starts, do_inter, rev_order, add_c4):
        lo, hi, d = span_delta(cur)
        blocks = []
        for s in starts:
            base = s * d - lo
            blocks.append([(l + base, r + base) for (l, r) in cur])
        S = []
        if do_inter:
            order = list(range(4))
            if rev_order: order.reverse()
            m = max(len(b) for b in blocks)
            for i in range(m):
                for idx in order:
                    blk = blocks[idx]
                    if i < len(blk):
                        S.append(blk[i])
        else:
            if rev_order:
                blocks = list(reversed(blocks))
            for blk in blocks:
                S.extend(blk)
        append_connectors(S, starts, d, add_cross4=add_c4)
        return S

    def build_micro(cur, budget, rid, alt):
        if budget <= 8 or not cur:
            return []
        glo, ghi, d0 = span_delta(cur)
        G = max(1, ghi - glo)
        # thin seed
        max_seed = 40 if alt else 32
        seed_sz = max(8, min(max_seed, len(cur) // 300))
        step = max(1, len(cur) // seed_sz)
        U = cur[::step][:seed_sz]
        if not U:
            return []
        ulo = min(l for l, r in U)
        # window fractions
        if not alt:
            window_fracs = [
                (0.12, 0.22), (0.35, 0.45),
                (0.58, 0.68), (0.80, 0.90),
            ]
        else:
            window_fracs = [
                (0.05, 0.15), (0.28, 0.38),
                (0.60, 0.70), (0.82, 0.92),
            ]
        # build blocks
        blocks = []
        for fa, fb in window_fracs:
            win_lo = glo + int(round(fa * G))
            base = win_lo - ulo
            b = [(l + base, r + base) for (l, r) in U]
            # diversify
            if alt and ((int(round(fa * 100)) // 5) % 2 == 0):
                b = list(reversed(b))
            blocks.append(b)
        # interleave
        micro = []
        m = max(len(b) for b in blocks)
        if alt:
            for i in range(m):
                for blk in reversed(blocks):
                    if i < len(blk):
                        micro.append(blk[i])
        else:
            for i in range(m):
                for blk in blocks:
                    if i < len(blk):
                        micro.append(blk[i])
        # local connectors at scale G
        append_connectors(micro, default_starts, max(1, G), add_cross4=(alt and cross4_enabled), lo_off=glo)
        # trim
        return micro[:budget]

    # Stage 1: KT spine
    for rid in range(int(rounds)):
        if 4 * len(T) + 4 > CAP:
            break
        starts = template_bank[rid % len(template_bank)] if rotate_starts else default_starts
        do_i = interleave_blocks and (rid % 2 == 0)
        rev = reverse_block_parity and (rid % 2 == 1)
        add4 = cross4_enabled and (rid == int(rounds) - 1)
        T = apply_round(T, starts, do_i, rev, add4)
        if len(T) >= CAP:
            return T[:CAP]

    if len(T) >= CAP - 8:
        return T[:CAP]

    # tail caps
    lo, hi, _ = span_delta(T)
    span = hi - lo
    def cap(a, b):
        L = lo + max(1, int(round(a * span)))
        R = lo + max(1, int(round(b * span)))
        if R <= L:
            R = L + 1
        return (L, R)
    for c in [cap(0.08, 0.60), cap(0.25, 0.75), cap(0.75, 0.92)]:
        if len(T) >= CAP:
            break
        T.append(c)

    # Stage 2: micro-phases
    steps = min(int(phase2_iters), 2)
    for i in range(steps):
        room = CAP - len(T)
        mi = build_micro(T, room, i, alt=False)
        if not mi:
            break
        T.extend(mi)
    # alternate micro
    room = CAP - len(T)
    if cross4_enabled and room > 8:
        mi = build_micro(T, room, steps, alt=True)
        if mi:
            T.extend(mi)

    return T[:CAP]

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()