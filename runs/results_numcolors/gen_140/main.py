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

    def micro_round(curr, rid, scale_div, extras):
        lo = min(l for l, _ in curr)
        hi = max(r for _, r in curr)
        G = max(1, hi - lo)
        delta = max(1, G // scale_div)
        # sample
        per = max(8, min(64, len(curr)//50))
        U = thin_seed(curr, per)
        if not U:
            return []
        # build blocks with parity‐reversal
        starts = (2, 6, 10, 14)
        ulo = min(l for l, _ in U)
        blocks = []
        for s in starts:
            base = s * delta - ulo
            blk = [(l+base, r+base) for (l, r) in U]
            if (s//2)%2 == rid%2:
                blk.reverse()
            blocks.append(blk)
        # interleave or reverse order
        out = []
        order = blocks if rid%2==0 else list(reversed(blocks))
        mlen = max(len(b) for b in blocks)
        for i in range(mlen):
            for b in order:
                if i < len(b):
                    out.append(b[i])
        # connectors: classic + extras
        s0, s1, s2, s3 = starts
        con = [
            ((s0-1)*delta+lo, (s1-1)*delta+lo),
            ((s2+2)*delta+lo, (s3+2)*delta+lo),
            ((s0+2)*delta+lo, (s2-1)*delta+lo),
            ((s1+2)*delta+lo, (s3-1)*delta+lo),
        ]
        # extras: list of (i, da, j, db) => (starts[i]+da, starts[j]+db)
        for (i, da, j, db) in extras:
            con.append(((starts[i]+da)*delta+lo, (starts[j]+db)*delta+lo))
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
        mr = micro_round(intervals, rid, p['scale_div'], p['extras'])
        if mr:
            intervals.extend(mr[:bud])

    # Phase 3: final sparse caps
    if len(intervals) < CAP:
        lo = min(l for l, _ in intervals)
        hi = max(r for _, r in intervals)
        G = max(1, hi - lo)
        d = max(1, G//3)
        caps = [
            (lo + d, lo + 4*d),
            (hi - 5*d, hi - 2*d),
            ((lo+hi)//2 - 2*d, (lo+hi)//2 + 3*d),
        ]
        rem = CAP - len(intervals)
        intervals.extend(caps[:rem])

    return intervals

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()