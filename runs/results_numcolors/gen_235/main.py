# EVOLVE-BLOCK-START

def construct_intervals(seed_count=1):
    """
    Same I/O as original: returns list of (l, r) open intervals.
    Implements an improved KT spine + multi-phase micro-gadget strategy.
    """
    # hard cap and gating margin
    CAP = 9800
    CAP_MARGIN = 32
    # density multiplier for spine
    K_DENSITY = 2

    # four strong KT templates
    template_bank = [
        (2, 6, 10, 14),
        (1, 5, 9, 13),
        (3, 7, 11, 15),
        (4, 8, 12, 16),
    ]

    # seed
    T = [(0,1)] if seed_count<=1 else [(i*3,i*3+1) for i in range(seed_count)]

    # helper: span and delta
    def span_delta(X):
        lo = min(l for l,r in X)
        hi = max(r for l,r in X)
        d = hi-lo if hi>lo else 1
        return lo, hi, d

    # append connectors
    def append_connectors(S, starts, delta, cross4=False):
        s0,s1,s2,s3 = starts
        S.append(((s0-1)*delta, (s1-1)*delta))
        S.append(((s2+2)*delta, (s3+2)*delta))
        S.append(((s0+2)*delta, (s2-1)*delta))
        S.append(((s1+2)*delta, (s3-1)*delta))
        if cross4:
            S.append(((s0+4)*delta, (s3+4)*delta))

    # Stage1: KT rounds with density K_DENSITY
    def apply_round(X, starts, interleave, rev, cross4=False):
        lo,hi,delta = span_delta(X)
        blocks=[]
        for s in starts:
            base = s*K_DENSITY*delta - lo
            blocks.append([(l+base, r+base) for l,r in X])
        S=[]
        if interleave:
            maxlen = max(len(b) for b in blocks)
            order = list(range(4))
            if rev: order.reverse()
            for i in range(maxlen):
                for idx in order:
                    blk=blocks[idx]
                    if i<len(blk): S.append(blk[i])
        else:
            if rev: blocks.reverse()
            for blk in blocks: S.extend(blk)
        append_connectors(S, starts, delta, cross4)
        return S

    # predict how many rounds fit
    def fit_rounds(initial, maxr):
        sz=initial; cnt=0
        for i in range(maxr):
            nxt = 4*sz+4
            if nxt>CAP: break
            sz=nxt; cnt+=1
        return cnt

    # run spine
    rounds=6
    depth = fit_rounds(len(T), rounds)
    for i in range(depth):
        starts = template_bank[i%4]
        inter = (i%2==0)
        rev   = (i%2==1)
        # sprinkle a long‐range connector on the last round
        cross4 = (i==depth-1)
        T = apply_round(T, starts, inter, rev, cross4=cross4)

    # if near capacity, skip micro-phases
    if len(T) > CAP - CAP_MARGIN:
        return T

    lo,hi,delta = span_delta(T)
    # simple long‐caps near tail
    def mkcap(a,b):
        L=lo+max(1,int(round(a*delta)))
        R=lo+max(1,int(round(b*delta)))
        return (L,R) if R>L else (L,L+1)
    tail_caps = [mkcap(0.08,0.60), mkcap(0.25,0.75), mkcap(0.75,0.92)]
    # insert near end to hit many active colors
    out = T[:]
    for idx,cap in enumerate(tail_caps):
        pos = len(out)-(2*idx+1)
        if pos<0: out.append(cap)
        else:     out.insert(pos, cap)
    T = out[:CAP]

    # Stage2: two delta2 micro-phases
    def build_micro(X, budget, windows, iter_id):
        glo,ghi,G = *span_delta(X), None
        glo,ghi,_ = span_delta(X)
        G = ghi-glo if ghi>glo else 1
        # thin seed
        seed_sz = max(8, min(40, len(X)//250))
        stride = max(1, len(X)//seed_sz)
        U = [X[j] for j in range(0,len(X),stride)][:seed_sz]
        if not U: return []
        ulo = min(l for l,r in U)
        blocks=[]
        for w_idx,(fa,fb) in enumerate(windows):
            win_lo = glo+int(round(fa*G))
            base = win_lo-ulo
            blk = [(l+base, r+base) for l,r in U]
            # break symmetry
            if (w_idx+iter_id)%2: blk.reverse()
            blocks.append(blk)
        # interleave
        micro=[]
        maxlen=max(len(b) for b in blocks)
        order=list(range(len(blocks)))
        if iter_id%2: order.reverse()
        for i in range(maxlen):
            for o in order:
                if i<len(blocks[o]): micro.append(blocks[o][i])
        # connectors
        con_fracs = [(0.08,0.30),(0.26,0.56),(0.44,0.78),(0.60,0.92)]
        for (a,b) in con_fracs:
            A=glo+int(round(a*G)); B=glo+int(round(b*G))
            if B>A: micro.append((A,B))
        return micro[:budget]

    # primary windows
    primary_w = [(0.12,0.22),(0.35,0.45),(0.58,0.68),(0.80,0.90)]
    # alternate windows
    alt_w     = [(0.05,0.15),(0.28,0.38),(0.60,0.70),(0.82,0.92)]
    # run two phases if room
    for phase,(wins,shifted) in enumerate(((primary_w,False),(alt_w,True))):
        if len(T) > CAP - CAP_MARGIN: break
        room = CAP - len(T)
        micro = build_micro(T, room, wins, phase)
        T.extend(micro)
    # parabolic micro-phase if still room
    if len(T) <= CAP - CAP_MARGIN:
        lo,hi,delta=span_delta(T)
        par = []
        for (a,b) in [(0.15,0.60),(0.25,0.80),(0.55,0.95)]:
            A=lo+int(round(a*delta)); B=lo+int(round(b*delta))
            if B>A: par.append((A,B))
        T.extend(par[:CAP-len(T)])

    return T[:CAP]

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()