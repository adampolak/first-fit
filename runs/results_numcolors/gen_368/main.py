# EVOLVE-BLOCK-START

def construct_intervals(enable_alt_microphase=True):
    """
    Construct an interval sequence to maximize FirstFit/omega.
    Uses a 7-round braided KT backbone with per-round densification,
    8 long connectors, two micro-phases, plus tower pins and quarter connectors.
    """
    CAP = 9800

    # 6-template start-pattern bank for KT rounds
    template_bank = [
        (2, 6, 10, 14),   # classic
        (1, 5, 9, 13),    # left-shift
        (3, 7, 11, 15),   # right-shift
        (4, 8, 12, 16),   # stretched-right
        (2, 5, 9, 14),    # mixed
        (3, 6, 10, 13),   # mixed2
    ]

    # Initial seed
    T = [(0, 1)]

    def _span(seq):
        lo = min(l for l, _ in seq)
        hi = max(r for _, r in seq)
        d = hi - lo
        return lo, hi, d if d > 0 else 1

    def _apply_braided_round(curr, starts, ridx):
        lo, hi, delta = _span(curr)
        # translate blocks
        blocks = []
        for s in starts:
            base = s * delta - lo
            blocks.append([(l + base, r + base) for l, r in curr])
        # braided interleaving: stripe widths rotate
        stripe = 1 + (ridx % 3)
        order = list(range(4))
        rot = ridx % 4
        order = order[rot:] + order[:rot]
        S = []
        idxs = [0]*4
        more = True
        while more:
            more = False
            for bi in order:
                b = blocks[bi]
                i0 = idxs[bi]
                if i0 < len(b):
                    more = True
                    j = min(i0+stripe, len(b))
                    S.extend(b[i0:j])
                    idxs[bi] = j
        # append 4 classic connectors
        s0, s1, s2, s3 = starts
        S.append(((s0-1)*delta, (s1-1)*delta))
        S.append(((s2+2)*delta, (s3+2)*delta))
        S.append(((s0+2)*delta, (s2-1)*delta))
        S.append(((s1+2)*delta, (s3-1)*delta))
        return S

    def _densify_center(curr, cap, ridx):
        """Inject shrunk copies of a thin seed from the central 20% span."""
        lo, hi, d = _span(curr)
        w_lo = lo + int(0.4*d)
        w_hi = lo + int(0.6*d)
        # sample up to 12 intervals uniformly
        n = len(curr)
        k = min(12, max(4, n//200))
        stride = max(1, n//k)
        U = [iv for iv in curr[::stride] if iv[0]>=w_lo and iv[1]<=w_hi][:k]
        eps = max(1, d//200)
        dens = []
        for l, r in U:
            if r-l > 2*eps:
                dens.append((l+eps, r-eps))
        # cap guard
        room = cap - len(curr)
        return dens[:room]

    # Stage 1: 7 braided KT rounds with densification
    for ridx in range(7):
        if 4*len(T)+4 > CAP:
            break
        starts = template_bank[ridx % len(template_bank)]
        T = _apply_braided_round(T, starts, ridx)
        if len(T) >= CAP:
            T = T[:CAP]; break
        # per-round densify
        dens = _densify_center(T, CAP, ridx)
        if dens:
            need = CAP - len(T)
            T.extend(dens[:need])
            if len(T) >= CAP:
                T = T[:CAP]; break

    # Long-range connectors (8) to force late colors
    lo, hi, d = _span(T)
    fracs = [
        (0.05,0.50),(0.10,0.62),(0.18,0.78),(0.25,0.85),
        (0.32,0.92),(0.40,0.60),(0.48,0.88),(0.55,0.95)
    ]
    for a,b in fracs:
        if len(T)>=CAP: break
        L = lo+int(a*d); R = lo+int(b*d)
        if R<=L: R=L+1
        T.append((L,R))

    if len(T)>=CAP-16:
        return T[:CAP]

    # Micro-phase A: four windows
    def _micro(curr, cap, windows, seed_factor, ridx, alt=False):
        lo, hi, d = _span(curr)
        room = cap - len(curr)
        if room<=8: return []
        # thin seed
        n = len(curr)
        k = min(50, max(8, n//seed_factor))
        stride = max(1, n//k)
        U = curr[::stride][:k]
        ulo = min(l for l,_ in U)
        blocks = []
        for (fa, fb) in windows:
            win = lo+int(fa*d)
            base = win-ulo
            blk = [(l+base,r+base) for l,r in U]
            if alt and (int(fa*100)%2==1):
                blk.reverse()
            blocks.append(blk)
        # interleave
        M = []
        ml = max(len(b) for b in blocks)
        for i in range(ml):
            for blk in (reversed(blocks) if alt else blocks):
                if i<len(blk):
                    M.append(blk[i])
        # connectors
        for (fa,fb) in [(0.08,0.30),(0.44,0.78),(0.26,0.56)]:
            A = lo+int(fa*d); B=lo+int(fb*d)
            if B>A: M.append((A,B))
        if alt:
            A=lo+int(0.18*d); B=lo+int(0.84*d)
            if B>A: M.append((A,B))
        return M[:room]

    # run micro-phase A
    winA = [(0.12,0.22),(0.35,0.45),(0.58,0.68),(0.80,0.90)]
    mA = _micro(T, CAP, winA, seed_factor=250, ridx=0, alt=False)
    if mA: T.extend(mA)
    # micro-phase B
    if enable_alt_microphase and len(T)<CAP-8:
        winB = [(0.05,0.15),(0.28,0.38),(0.60,0.70),(0.82,0.92)]
        mB = _micro(T, CAP, winB, seed_factor=220, ridx=1, alt=True)
        if mB: T.extend(mB)

    # Tower pins: four interior windows, short eps-length intervals
    if len(T)<CAP-8:
        lo, hi, d = _span(T)
        eps = max(1, d//4096)
        wins = [0.25,0.45,0.65,0.85]
        pins = []
        for wi,w in enumerate(wins):
            base = lo+int(w*d)
            cnt = min((CAP-len(T))//4, 6)
            for j in range(cnt):
                L = base + j*(eps+2) + wi%3
                R = L+eps
                if R>L:
                    pins.append((L,R))
        need = CAP-len(T)
        T.extend(pins[:need])

    # Final quarter-span connectors
    if len(T)<CAP:
        lo, hi, d = _span(T)
        qs = [(0.25,0.75),(0.50-0.1,0.50+0.1),(0.75,0.95)]
        for (fa,fb) in qs:
            if len(T)>=CAP: break
            A=lo+int(fa*d); B=lo+int(fb*d)
            if B<=A: B=A+1
            T.append((A,B))

    return T[:CAP]

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()