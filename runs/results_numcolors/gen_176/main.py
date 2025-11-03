# EVOLVE-BLOCK-START

def construct_intervals(seed_count=1):
    # Hard capacity
    CAP = 9800

    # Extended start‐pattern bank
    start_patterns = [
        (2,6,10,14),(1,5,9,13),(3,7,11,15),(4,8,12,16),
        (2,4,8,12),(3,6,9,12),(1,4,7,10),(0,5,10,15)
    ]

    # Initial seed
    if seed_count <= 1:
        T = [(0,1)]
    else:
        # sparse multi‐seed up to 4 intervals
        step = max(1, seed_count)
        T = [(i*step, i*step+1) for i in range(seed_count)][:4]

    # Helpers
    def span_delta(U):
        lo = min(l for l,_ in U)
        hi = max(r for _,r in U)
        d = hi - lo
        return lo, hi, d if d>0 else 1

    def append_connectors(S, starts, delta, extra_cross4=False):
        s0,s1,s2,s3 = starts
        S.append(((s0-1)*delta, (s1-1)*delta))
        S.append(((s2+2)*delta, (s3+2)*delta))
        S.append(((s0+2)*delta, (s2-1)*delta))
        S.append(((s1+2)*delta, (s3-1)*delta))
        if extra_cross4:
            S.append(((s0+4)*delta, (s3+4)*delta))

    def feasible_depth(n, maxr):
        sz, d = n,0
        for _ in range(maxr):
            nxt = 4*sz + 4
            if nxt > CAP: break
            sz, d = nxt, d+1
        return d

    # Stage 1: Spine
    rounds = 6
    depth = feasible_depth(len(T), rounds)
    for ridx in range(depth):
        starts = start_patterns[ridx % len(start_patterns)]
        lo,hi,delta = span_delta(T)
        K = 2 if ridx >= depth-2 else 1
        # build blocks
        blocks = []
        for s in starts:
            base = s * delta * K - lo
            blocks.append([(l+base, r+base) for (l,r) in T])
        # interleaving policy
        inter = (ridx % 2 == 0)
        rev   = (ridx % 2 == 1)
        S = []
        if inter:
            order = list(range(len(blocks)))
            if rev: order.reverse()
            m = max(len(b) for b in blocks)
            for i in range(m):
                for idx in order:
                    if i < len(blocks[idx]):
                        S.append(blocks[idx][i])
        else:
            order = list(range(len(blocks)))
            if rev: order.reverse()
            for idx in order:
                S.extend(blocks[idx])
        # connectors
        extra = (ridx >= depth-2)
        append_connectors(S, starts, delta, extra_cross4=extra)
        T = S
        if len(T) >= CAP-16:
            break

    # Stage 2: Two delta₂ micro‐phases
    def build_micro(U, window_fracs, budget, iter_id):
        if budget <= 4 or not U: return []
        glo,ghi,G = *span_delta(U),
        # thin seed ≤32
        seed_sz = min(32, max(8, len(U)//300))
        step    = max(1, len(U)//seed_sz)
        V = [U[i] for i in range(0,len(U),step)][:seed_sz]
        if not V: return []
        ulo = min(l for l,_ in V)
        # blocks per window
        blocks = []
        for (fa,fb) in window_fracs:
            wlo = glo + int(round(fa*G))
            base = wlo - ulo
            blocks.append([(l+base, r+base) for (l,r) in V])
        # interleave or reverse‐interleave
        micro=[]
        inter = (iter_id % 2 == 0)
        order = list(range(4))
        if not inter: order.reverse()
        m = max(len(b) for b in blocks)
        for i in range(m):
            for idx in order:
                if i < len(blocks[idx]):
                    micro.append(blocks[idx][i])
        # connectors & caps
        # 4 classic + 1 cross4 + 3 safe caps
        s0,s1,s2,s3 = (2,6,10,14)
        d2 = max(1, G//2)
        # classic at delta2+glo
        micro.append(((s0-1)*d2+glo, (s1-1)*d2+glo))
        micro.append(((s2+2)*d2+glo, (s3+2)*d2+glo))
        micro.append(((s0+2)*d2+glo, (s2-1)*d2+glo))
        micro.append(((s1+2)*d2+glo, (s3-1)*d2+glo))
        micro.append(((s0+4)*d2+glo, (s3+4)*d2+glo))
        # three safe caps
        mid = glo + G//2
        caps = [
            (glo + int(0.08*G), glo + int(0.30*G)),
            (glo + int(0.60*G), glo + int(0.92*G)),
            (mid - max(1,d2//8), mid + max(1,d2//8))
        ]
        for c in caps:
            if c[1]>c[0]: micro.append(c)
        return micro[:budget]

    rem = CAP - len(T)
    if rem > 16:
        wf1 = [(0.10,0.20),(0.32,0.42),(0.55,0.65),(0.78,0.88)]
        m1 = build_micro(T, wf1, rem//2, 0)
        T.extend(m1[:CAP-len(T)])
    rem = CAP - len(T)
    if rem > 8:
        wf2 = [(0.05,0.15),(0.25,0.35),(0.50,0.60),(0.75,0.85)]
        m2 = build_micro(T, wf2, rem, 1)
        T.extend(m2[:CAP-len(T)])

    # final trim
    return T[:CAP]

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()