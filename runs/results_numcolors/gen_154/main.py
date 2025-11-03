# EVOLVE-BLOCK-START

def construct_intervals(seed_count=1):
    """
    Enhanced multi-scale spine + three-round micro-phase.
    Returns a list of (l,r) integer intervals for FirstFit.
    """
    MAX_CAP = 9800
    # Four 4-start templates
    spine_templates = [
        (2, 6, 10, 14),
        (1, 5, 9, 13),
        (3, 7, 11, 15),
        (4, 8, 12, 16),
    ]
    # Seed
    T = [(0.0, 1.0)]
    # Stage 1: 5 spine rounds
    for ridx in range(5):
        lo = min(l for l, r in T)
        hi = max(r for l, r in T)
        delta = hi - lo if hi > lo else 1.0
        starts = spine_templates[ridx % len(spine_templates)]
        K = 1 + (ridx % 2)  # alternate densification
        # build 4 blocks
        blocks = []
        for s in starts:
            base = (s * K) * delta - lo
            blocks.append([(l + base, r + base) for (l, r) in T])
        # parity interleaving: even rounds interleave, odd sequential
        S = []
        if ridx % 2 == 0:
            maxlen = max(len(b) for b in blocks)
            for i in range(maxlen):
                for blk in blocks:
                    if i < len(blk):
                        S.append(blk[i])
        else:
            for blk in blocks:
                S.extend(blk)
        # classic 4 connectors
        s0, s1, s2, s3 = starts
        for a, b in [
            ((s0 - 1) * K * delta, (s1 - 1) * K * delta),
            ((s2 + 2) * K * delta, (s3 + 2) * K * delta),
            ((s0 + 2) * K * delta, (s2 - 1) * K * delta),
            ((s1 + 2) * K * delta, (s3 - 1) * K * delta),
        ]:
            S.append((a, b))
        T = S
        if len(T) >= MAX_CAP - 16:
            break

    # If spine already near cap, convert and return
    if len(T) >= MAX_CAP - 100:
        # simple int conversion
        out = []
        for (l, r) in T:
            li = int(round(l))
            ri = int(round(r))
            if ri <= li: ri = li + 1
            out.append((li, ri))
            if len(out) >= 10000:
                break
        return out

    # Stage 2: three micro rounds
    def thin_seed(U, m):
        if not U or m <= 0:
            return []
        step = max(1, len(U) // m)
        return U[::step][:m]

    def micro_round(curr, rid, budget):
        if budget <= 0 or not curr:
            return []
        glo = min(l for l, r in curr)
        ghi = max(r for l, r in curr)
        G = max(1.0, ghi - glo)
        delta2 = G / 2.0
        # sample
        per = max(8, min(96, budget // 8))
        U = thin_seed(curr, per)
        if not U:
            return []
        # template & density
        starts = spine_templates[rid % len(spine_templates)]
        K2 = 1 + (rid % 2)
        blocks = []
        # preload blockers params
        pb = 2
        bl = 0.05 * G
        bg = 0.03 * G
        for idx, s in enumerate(starts):
            base = (s * K2) * delta2 - min(l for l, r in U)
            block = []
            # small preload blockers
            for b in range(pb):
                st = base + (0.1 + b * (bl + bg))
                block.append((glo + st, glo + st + bl))
            # clone (with optional reversal)
            seq = reversed(U) if (idx % 2 == rid % 2) else U
            for (l, r) in seq:
                block.append((l + base, r + base))
            blocks.append(block)
        # interleave or sequential
        M = []
        if rid % 2 == 0:
            maxlen = max(len(b) for b in blocks)
            for i in range(maxlen):
                for blk in blocks:
                    if i < len(blk):
                        M.append(blk[i])
        else:
            maxlen = max(len(b) for b in blocks)
            for i in range(maxlen):
                for blk in reversed(blocks):
                    if i < len(blk):
                        M.append(blk[i])
        # micro-bridges
        for i in range(len(starts) - 1):
            M.append((glo + (starts[i] + 0.6) * delta2,
                      glo + (starts[i + 1] - 0.6) * delta2))
        # connectors cross1..cross4 + cross3
        s0, s1, s2, s3 = starts
        for a, b in [
            ((s0 - 1) * K2 * delta2 + glo, (s1 - 1) * K2 * delta2 + glo),
            ((s2 + 2) * K2 * delta2 + glo, (s3 + 2) * K2 * delta2 + glo),
            ((s0 + 2) * K2 * delta2 + glo, (s2 - 1) * K2 * delta2 + glo),
            ((s1 + 2) * K2 * delta2 + glo, (s3 - 1) * K2 * delta2 + glo),
            ((s0 + 3) * K2 * delta2 + glo, (s3 + 3) * K2 * delta2 + glo),
            ((s1 + 4) * K2 * delta2 + glo, (s2 + 4) * K2 * delta2 + glo),
        ]:
            M.append((a, b))
        return M[:budget]

    remaining = MAX_CAP - len(T)
    for rid in range(3):
        bud = remaining // (3 - rid)
        mr = micro_round(T, rid, bud)
        room = MAX_CAP - len(T)
        if mr:
            T.extend(mr[:room])
        remaining = MAX_CAP - len(T)

    # Final conversion to integer intervals
    if not T:
        return []
    minl = min(l for l, r in T)
    sf = 1000.0
    out = []
    for (l, r) in T:
        li = int(round((l - minl) * sf))
        ri = int(round((r - minl) * sf))
        if ri <= li:
            ri = li + 1
        out.append((li, ri))
        if len(out) >= 10000:
            break
    return out

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()