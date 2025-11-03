# EVOLVE-BLOCK-START

def construct_intervals():
    """
    Adaptive KT backbone with per‐round densification and two‐phase micro‐budgeting.
    Returns a list of open intervals (l, r).
    """
    CAP = 9800

    # 6‐template backbone for richer mixing
    template_bank = [
        (2, 6, 10, 14),
        (1, 5, 9, 13),
        (3, 7, 11, 15),
        (4, 8, 12, 16),
        (5, 9, 13, 17),
        (0, 7, 14, 21),
    ]

    # seed
    T = [(0, 1)]

    def span_info(S):
        lo = min(l for l, r in S)
        hi = max(r for l, r in S)
        d = hi - lo
        return lo, hi, max(d, 1)

    def apply_round(S, starts, interleave, reverse):
        lo, hi, d = span_info(S)
        blocks = []
        for s in starts:
            base = s * d - lo
            blocks.append([(l + base, r + base) for l, r in S])
        R = []
        if interleave:
            order = list(range(4))
            if reverse:
                order.reverse()
            m = max(len(b) for b in blocks)
            for i in range(m):
                for idx in order:
                    blk = blocks[idx]
                    if i < len(blk):
                        R.append(blk[i])
        else:
            seq = blocks[::-1] if reverse else blocks
            for blk in seq:
                R.extend(blk)
        # classic connectors
        s0, s1, s2, s3 = starts
        R.append(((s0 - 1) * d, (s1 - 1) * d))
        R.append(((s2 + 2) * d, (s3 + 2) * d))
        R.append(((s0 + 2) * d, (s2 - 1) * d))
        R.append(((s1 + 2) * d, (s3 - 1) * d))
        return R

    # Stage 1: up to 6 rounds (4*S+4 growth) under CAP
    for ridx in range(6):
        if 4 * len(T) + 4 > CAP:
            break
        starts = template_bank[ridx % len(template_bank)]
        inter = (ridx % 2 == 0)
        rev = (ridx % 3 == 2)
        T = apply_round(T, starts, interleave=inter, reverse=rev)

        # CAP‐aware per‐round densification: shrink & reinsert small sample
        if len(T) < CAP - 5:
            lo, hi, d = span_info(T)
            eps = max(1, d // 800)
            sample_sz = min(5, max(1, len(T) // 2000))
            for l, r in T[:sample_sz]:
                if r - l > 2 * eps:
                    T.append((l + eps, r - eps))
                    if len(T) >= CAP:
                        break

    # first micro‐budget
    room = CAP - len(T)
    if room <= 8:
        return T[:CAP]
    b1 = int(room * 0.6)
    b2 = room - b1

    def build_micro(S, budget, windows):
        lo, hi, G = span_info(S)
        seed_sz = min(16, max(8, len(S) // 200))
        stride = max(1, len(S) // seed_sz)
        U = [S[i] for i in range(0, len(S), stride)][:seed_sz]
        out = []
        # build blocks
        for fa, fb in windows:
            win_lo = lo + int(round(fa * G))
            base = win_lo - min(l for l, r in U)
            blk = [(l + base, r + base) for l, r in U]
            out.extend(blk)
        # connectors
        for fa, fb in windows[:4]:
            a = lo + int(round(fa * G - 0.02 * G))
            b = lo + int(round(fb * G + 0.02 * G))
            if b > a:
                out.append((a, b))
        return out[:budget]

    # Micro-phase A: 4-window
    winA = [(0.12,0.22),(0.36,0.46),(0.60,0.70),(0.84,0.94)]
    micro1 = build_micro(T, b1, winA)
    T.extend(micro1)

    # deterministic long‐range connectors
    lo, hi, d = span_info(T)
    long_conns = [
        (lo + int(round(0.05 * d)), hi - int(round(0.05 * d))),
        (lo + int(round(0.20 * d)), lo + int(round(0.85 * d))),
        (lo + int(round(0.10 * d)), hi - int(round(0.20 * d))),
    ]
    for a, b in long_conns:
        if len(T) < CAP and b > a:
            T.append((a, b))

    # Micro-phase B: 5-window alternate
    room2 = CAP - len(T)
    if room2 > 8 and b2 > 0:
        winB = [(0.06,0.16),(0.28,0.38),(0.50,0.60),(0.72,0.82),(0.88,0.94)]
        micro2 = build_micro(T, min(b2, room2), winB)
        T.extend(micro2)

    return T[:CAP]

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()