# EVOLVE-BLOCK-START

def construct_intervals():
    """
    Construct a sequence of open intervals for FirstFit,
    maximizing FF colors / clique number via an improved multi-phase strategy.
    Returns:
      intervals: list of (l, r) tuples.
    """
    CAP = 9800

    # Phase 1: KT‐style spine with rotated templates and alternating multiplier K
    template_bank = [
        (2, 6, 10, 14),  # A: classic
        (1, 5, 9, 13),   # B: left-shifted
        (3, 7, 11, 15),  # C: right-shifted
        (4, 8, 12, 16),  # D: alternate
    ]
    K_values = [1, 2]
    spine_rounds = 6

    intervals = [(0, 1)]  # seed

    def apply_spine_round(curr, starts, K):
        lo = min(l for l, _ in curr)
        hi = max(r for _, r in curr)
        delta = max(1, hi - lo)
        out = []
        # four translated blocks
        for s in starts:
            base = s * delta * K - lo
            for (l, r) in curr:
                out.append((l + base, r + base))
        # classic connectors
        s0, s1, s2, s3 = starts
        con = [
            ((s0-1)*delta*K, (s1-1)*delta*K),
            ((s2+2)*delta*K, (s3+2)*delta*K),
            ((s0+2)*delta*K, (s2-1)*delta*K),
            ((s1+2)*delta*K, (s3-1)*delta*K),
        ]
        out.extend(con)
        return out

    # build spine
    for i in range(spine_rounds):
        if 4*len(intervals)+4 > CAP:
            break
        tpl = template_bank[i % len(template_bank)]
        K = K_values[i % len(K_values)]
        intervals = apply_spine_round(intervals, tpl, K)
    if len(intervals) >= CAP - 16:
        return intervals

    # Phase 2: three micro‐rounds at scale divisors 2,3,4
    def thin_seed(curr, max_seed):
        n = len(curr)
        if n == 0 or max_seed <= 0:
            return []
        step = max(1, n // max_seed)
        return curr[::step][:max_seed]

    def micro_round(curr, rid, scale_div, extras, budget_hint):
        lo = min(l for l, _ in curr)
        hi = max(r for _, r in curr)
        G = max(1, hi - lo)
        delta = max(1, G // scale_div)

        # adaptive sample size: slightly larger sample to strengthen micro pressure
        per_floor = 16
        per = min(192, max(per_floor, budget_hint // 10, len(curr) // 64))
        U = thin_seed(curr, per)
        if not U:
            return []

        # build blocks with parity-reversal and deterministic rotation; then shrink
        starts = (2, 6, 10, 14)
        ulo = min(l for l, _ in U)
        eps = max(1, delta // 64)  # micro pin length
        blocks = []
        for idx, s in enumerate(starts):
            base = s * delta - ulo
            raw_blk = [(l + base, r + base) for (l, r) in U]
            if (s // 2) % 2 == rid % 2:
                raw_blk.reverse()

            # two localized blockers to preoccupy small colors inside the block
            blockers = []
            if raw_blk:
                left_anchor = int(raw_blk[0][0] + max(1, delta // 10))
                right_anchor = int(raw_blk[-1][1] - max(1, delta // 8))
                blockers.append((left_anchor, left_anchor + eps))
                if right_anchor - eps > left_anchor:
                    blockers.append((right_anchor, right_anchor + eps))

            # shrink each interval in the block to a short "pin" to limit omega
            shrunk = []
            for (l, r) in raw_blk:
                mid = (l + r) // 2
                L = mid - eps // 2
                R = L + eps
                if R > L:
                    shrunk.append((L, R))

            # present blockers first in order to consume smaller colors early
            if blockers:
                shrunk = blockers + shrunk
            blocks.append(shrunk)

        # Interleave with a rotated order to improve cross-block coupling
        out = []
        order_idx = list(range(len(blocks)))
        krot = rid % len(order_idx)
        order_idx = order_idx[krot:] + order_idx[:krot]
        mlen = max(len(b) for b in blocks)
        for i in range(mlen):
            for idx in order_idx:
                b = blocks[idx]
                if i < len(b):
                    out.append(b[i])

        # connectors: keep only the two cap connectors at micro scale to avoid omega blow-up
        s0, s1, s2, s3 = starts
        con = [
            ((s0 - 1) * delta + lo, (s1 - 1) * delta + lo),  # left cap
            ((s2 + 2) * delta + lo, (s3 + 2) * delta + lo),  # right cap
        ]
        out.extend(con)
        return out

    subphases = [
        {'scale_div':2, 'extras':[(0,3,3,3)]},  # cross3
        {'scale_div':3, 'extras':[(1,4,2,4)]},  # cross4
        {'scale_div':4, 'extras':[(0,5,3,5)]},  # cross5
    ]

    for rid, p in enumerate(subphases):
        if len(intervals) >= CAP:
            break
        bud = CAP - len(intervals)
        # pass remaining-budget hint so micro_round can adapt sampling density and extras
        mr = micro_round(intervals, rid, p['scale_div'], p['extras'], budget_hint=bud)
        if mr:
            intervals.extend(mr[:bud])

    # Phase 3: final micro-pins (short caps) to avoid inflating omega
    if len(intervals) < CAP:
        lo = min(l for l, _ in intervals)
        hi = max(r for _, r in intervals)
        G = max(1, hi - lo)
        step = max(1, G // 9)
        eps = max(1, step // 8)
        centers = [lo + 2 * step, (lo + hi) // 2, hi - 2 * step, lo + 5 * step]
        pins = []
        for c in centers:
            L = int(c - eps // 2)
            R = L + int(eps)
            if R > L:
                pins.append((L, R))
        rem = CAP - len(intervals)
        intervals.extend(pins[:rem])

    return intervals

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()