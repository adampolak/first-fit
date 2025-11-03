# EVOLVE-BLOCK-START

from math import gcd
# Basic overlap, FirstFit and clique routines
def overlaps(a,b):
    (l1,r1),(l2,r2)=a,b
    return max(l1,l2)<min(r1,r2)

def firstfit_colors(intervals):
    last_end=[]
    for l,r in intervals:
        placed=False
        for i,le in enumerate(last_end):
            if l>=le:
                last_end[i]=r
                placed=True
                break
        if not placed:
            last_end.append(r)
    return len(last_end)

def clique_number(intervals):
    events=[]
    for l,r in intervals:
        if l<r:
            events.append((l,1)); events.append((r,-1))
    events.sort(key=lambda e:(e[0],e[1]))
    cur=best=0
    for _,d in events:
        cur+=d
        if cur>best: best=cur
    return best

# Two‐stage pruning
def prune_stage1(T, tc, to):
    changed=True
    while changed:
        changed=False
        # try removing longest intervals first
        for i,iv in sorted(enumerate(T), key=lambda x:-(x[1][1]-x[1][0])):
            U=T[:i]+T[i+1:]
            if firstfit_colors(U)==tc and clique_number(U)==to:
                T=U; changed=True; break
    return T

def prune_stage2(T, tc, to):
    # find a witness point for ω
    events=[]
    for l,r in T:
        events.append((l,1)); events.append((r,-1))
    events.sort(key=lambda e:(e[0],e[1]))
    cur=best=0; wp=events[0][0]
    for x,d in events:
        cur+=d
        if cur>best:
            best=cur; wp=x
    cover={i for i,(l,r) in enumerate(T) if l<wp<r}
    changed=True
    while changed:
        changed=False
        for i,iv in sorted(enumerate(T), key=lambda x:-(x[1][1]-x[1][0])):
            if i in cover: continue
            U=T[:i]+T[i+1:]
            if firstfit_colors(U)==tc and clique_number(U)==to:
                T=U; changed=True; break
    return T

def normalize(T):
    if not T: return []
    pts=sorted({x for seg in T for x in seg})
    mp={v:i*2 for i,v in enumerate(pts)}
    L=[(mp[l],mp[r]) for l,r in T]
    mn=min(l for l,_ in L)
    L=[(l-mn,r-mn) for l,r in L]
    g=0
    for a,b in L:
        g=gcd(g,abs(a)); g=gcd(g,abs(b))
    if g>1:
        L=[(a//g,b//g) for a,b in L]
    return L

# Wave templates
SHORT_WAVE=[(1,3),(3,5),(5,7)]
LONG_WAVE=[(8,12),(12,16),(16,20)]
# Four tiling patterns
TILING=[(2,6,10,14),(1,5,9,13),(3,7,11,15),(0,4,8,12)]
# Memoization cache
_eval_cache={}

def evaluate_cached(intervals):
    key=tuple(normalize(intervals))
    if key in _eval_cache:
        return _eval_cache[key]
    Tn=normalize(intervals)
    om=clique_number(Tn)
    cols=firstfit_colors(Tn)
    score=cols/om - 1e-6*(len(Tn)/10000.0)
    _eval_cache[key]=(score,om,cols,len(Tn),Tn)
    return _eval_cache[key]

def build_pattern(base, depth, blockers, translation, anchor, wave_type, eps):
    T=list(base)
    for lvl in range(depth):
        lo=min(l for l,r in T); hi=max(r for l,r in T)
        d=hi-lo; c=(lo+hi)/2.0
        # cycle tiling with perturbation
        offs=[s+eps for s in TILING[lvl%4]]
        S=[]
        for s in offs:
            off = (d*s - lo) if translation=='left' else (d*s - c)
            for l,r in T:
                S.append((l+off,r+off))
        # blockers
        for a,b in blockers:
            if anchor=='left':
                S.append((d*(a+eps),d*(b+eps)))
            else:
                S.append((d*(a+eps)-c,d*(b+eps)-c))
        # wave augmentation every other level
        wave = SHORT_WAVE if wave_type=='short' else LONG_WAVE
        if lvl%2==0:
            for a,b in wave:
                S.append((d*a,d*b))
        T=S
    return T

def construct_intervals():
    base_seeds=[[(0.0,1.0)],[(0.0,1.0),(2.0,3.0)]]
    blockers_set=[((1,5),(12,16),(4,9),(8,13)), ((0,4),(6,10),(8,12),(14,18))]
    translations=['left','center']
    anchors=['left','center']
    wave_types=['short','long']
    epsilons=[-0.1,0.0,0.1]
    depths=[3,4,5,6]

    best=None
    for base in base_seeds:
        for depth in depths:
            # cap instance size
            if (4**depth)*(len(base)+len(SHORT_WAVE))>3000: continue
            for blockers in blockers_set:
                for translation in translations:
                    for anchor in anchors:
                        for wave_type in wave_types:
                            for eps in epsilons:
                                T=build_pattern(base,depth,blockers,translation,anchor,wave_type,eps)
                                sc,om,cols,n,_=evaluate_cached(T)
                                T1=prune_stage1(T,cols,om)
                                T2=prune_stage2(T1,cols,om)
                                cand=evaluate_cached(T2)
                                if best is None or \
                                   cand[0]>best[0]+1e-9 or \
                                   (abs(cand[0]-best[0])<1e-9 and (cand[3]<best[3] or (cand[3]==best[3] and cand[2]>best[2]))):
                                    best=cand
    if best:
        return best[4]
    # fallback classic
    T=[(0.0,1.0)]
    for _ in range(4):
        lo=min(l for l,r in T); hi=max(r for l,r in T); d=hi-lo
        S=[(l+d*s-lo,r+d*s-lo) for s in (2,6,10,14) for l,r in T]
        S+=[(d*1,d*5),(d*12,d*16),(d*4,d*9),(d*8,d*13)]
        T=S
    return normalize(T)

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()