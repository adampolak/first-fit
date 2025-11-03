# EVOLVE-BLOCK-START

from math import gcd

# --- Normalization and Utilities ---

def _normalize_grid(intervals):
    """
    Normalize endpoints to a compact integer grid while preserving order.
    Each unique endpoint is mapped to an increasing even integer.
    """
    endpoints = sorted(set(x for seg in intervals for x in seg))
    coord = {}
    cur = 0
    for e in endpoints:
        coord[e] = cur
        cur += 2
    return [(coord[l], coord[r]) for l, r in intervals]

def _ff_count(intervals):
    """
    Fast FirstFit color count using per-color last end tracking.
    """
    last_end = []
    for l, r in intervals:
        placed = False
        for i, le in enumerate(last_end):
            if l >= le:
                last_end[i] = r
                placed = True
                break
        if not placed:
            last_end.append(r)
    return len(last_end)

def _omega_open(intervals):
    """
    Compute omega (max intervals covering a point) via sweep.
    Process right endpoints before left for open intervals.
    """
    events = []
    for l, r in intervals:
        if l < r:
            events.append((l, 1))
            events.append((r, -1))
    events.sort(key=lambda e:(e[0], 0 if e[1]==-1 else 1))
    cur = best = 0
    for _, d in events:
        cur += d
        if cur > best:
            best = cur
    return best

def overlaps(a, b):
    """Test open-interval overlap."""
    return max(a[0], b[0]) < min(a[1], b[1])

# --- Baseline Generator (_build_candidate) from inspiration ---

def _build_candidate(depth, starts, blockers, schedule='after', extra_first=False, extra_last=False, translation='left'):
    """
    Build recursive baseline with scaled copies and blockers.
    """
    T = [(0.0, 1.0)]
    for i in range(depth):
        lo = min(l for l, r in T)
        hi = max(r for l, r in T)
        delta = hi - lo
        center = (lo + hi)/2.0
        offs = tuple(starts)
        if extra_first and i==0 and 18 not in offs:
            offs = offs + (18,)
        if extra_last and i==depth-1 and 18 not in offs:
            offs = offs + (18,)

        def make_copies(src, offs):
            S=[]
            for st in offs:
                offset = delta*st - (lo if translation=='left' else center)
                for l,r in src:
                    S.append((l+offset, r+offset))
            return S

        blk = []
        for a,b in blockers:
            o = delta*a - (0 if translation=='left' else center)
            p = delta*b - (0 if translation=='left' else center)
            blk.append((o,p))

        if schedule=='after':
            S = make_copies(T,offs) + blk
        elif schedule=='before':
            S = blk + make_copies(T,offs)
        elif schedule=='split':
            h = len(offs)//2
            S = make_copies(T,offs[:h]) + blk + make_copies(T,offs[h:])
        elif schedule=='interleaved':
            S=[]
            for idx,st in enumerate(offs):
                offset = delta*st - (lo if translation=='left' else center)
                for l,r in T:
                    S.append((l+offset,r+offset))
                if idx<len(blk):
                    S.append(blk[idx])
        else:
            S = make_copies(T,offs) + blk

        T=S
    return T

# --- Augmentation: corridor waves from current implementation ---

def firstfit_partition(intervals):
    """
    Partition intervals by FirstFit using overlaps.
    """
    colors=[]
    for iv in intervals:
        placed=False
        for cls in colors:
            if not any(overlaps(iv,u) for u in cls):
                cls.append(iv)
                placed=True
                break
        if not placed:
            colors.append([iv])
    return colors

def coverage_cells(intervals):
    """Return cells (l,r,count) for coverage between endpoints."""
    events=[]
    for l,r in intervals:
        if l<r:
            events.append((l,1))
            events.append((r,-1))
    if not events:
        return []
    events.sort(key=lambda e:(e[0],0 if e[1]==-1 else 1))
    xs = sorted({x for x,_ in events})
    idx=0
    cur=0
    counts={}
    for x in xs:
        while idx<len(events) and events[idx][0]==x:
            cur+=events[idx][1]
            idx+=1
        counts[x]=cur
    cells=[]
    for i in range(len(xs)-1):
        a,b=xs[i],xs[i+1]
        m=counts[a]
        if a<b:
            cells.append((a,b,m))
    return cells

def extract_runs(cells, cap):
    """
    Extract runs where coverage <= cap.
    """
    runs=[]
    curL=curR=None
    cur_min=1e18
    for l,r,m in cells:
        if m<=cap:
            if curL is None:
                curL,curR,cur_min = l,r,m
            elif abs(l-curR)<1e-12:
                curR=r
                cur_min=min(cur_min,m)
            else:
                runs.append((curL,curR,cur_min))
                curL,curR,cur_min = l,r,m
        else:
            if curL is not None:
                runs.append((curL,curR,cur_min))
                curL=None
    if curL is not None:
        runs.append((curL,curR,cur_min))
    return runs

def colors_hit_by_segment(seg, color_classes):
    """Count how many color classes are hit by segment."""
    L,R=seg
    cnt=0
    for cls in color_classes:
        if any(max(l,L)<min(r,R) for l,r in cls):
            cnt+=1
    return cnt

def pick_augmentable_run(baseline, omega_cap):
    """
    Pick a run [L,R] covering all colors and with room for two waves.
    """
    cells=coverage_cells(baseline)
    runs=extract_runs(cells, omega_cap-1)
    classes=firstfit_partition(baseline)
    total=len(classes)
    best=None
    for L,R,min_m in runs:
        if R<=L: continue
        hit=colors_hit_by_segment((L,R),classes)
        if hit<total: continue
        allow_two = (min_m <= omega_cap-2)
        core=None
        if allow_two:
            best_len=0
            for a,b,m in cells:
                if m<=omega_cap-2:
                    lo,hi=max(L,a),min(R,b)
                    if hi>lo and hi-lo>best_len:
                        best_len=hi-lo
                        core=(lo+1e-3,hi-1e-3)
            allow_two = core is not None
        if not allow_two:
            mid=(L+R)/2
            core=(mid-1e-3,mid+1e-3)
        score=(1 if allow_two else 0, R-L)
        if best is None or score>best[0]:
            best=(score,(L,R,allow_two,core))
    return best[1] if best else None

def add_wave(intervals, seg, omega_cap):
    """
    Try adding seg; require new color and omega<=cap.
    """
    cls=firstfit_partition(intervals)
    if colors_hit_by_segment(seg,cls)<len(cls):
        return intervals,False
    cand=intervals+[seg]
    if _omega_open(cand)>omega_cap:
        return intervals,False
    if _ff_count(cand)<=_ff_count(intervals):
        return intervals,False
    return cand,True

def augment_with_corridor_waves(base):
    """
    Add up to two waves to base.
    """
    if not base: return base
    omega_cap=_omega_open(base)
    if omega_cap<=1: return base
    pick=pick_augmentable_run(base,omega_cap)
    if not pick: return base
    L,R,allow_two,core=pick
    c_lo,c_hi=core
    c_lo=max(L+1e-3,min(c_lo,R-1e-3))
    c_hi=min(R-1e-3,max(c_hi,L+1e-3))
    if not (L<c_lo<c_hi<R): return base
    J1=(L,c_hi)
    T1,ok1=add_wave(base,J1,omega_cap)
    if not ok1:
        J1=(c_lo,R)
        T1,ok1=add_wave(base,J1,omega_cap)
        if not ok1:
            return base
    J2=(c_lo,R) if J1[0]==L else (L,c_hi)
    T2,ok2=add_wave(T1,J2,omega_cap)
    if not ok2 and allow_two:
        for f in [0.1,0.2,0.3]:
            if J1[0]==L:
                J2_try=(c_lo+f*(c_hi-c_lo),R)
            else:
                J2_try=(L,c_hi-f*(c_hi-c_lo))
            T2,ok2=add_wave(T1,J2_try,omega_cap)
            if ok2: break
    return T2 if ok2 else T1

# --- Orchestrator ---

def construct_intervals():
    """
    Build baseline via _build_candidate then augment.
    """
    offsets_set=[(2,6,10,14),(1,5,9,13),(3,7,11,15),(0,4,8,12)]
    blockers=[(1,5),(12,16),(4,9),(8,13)]
    k=4
    best=None
    for offs in offsets_set:
        base=_build_candidate(k, offs, blockers, schedule='after',
                              extra_first=False, extra_last=False, translation='left')
        norm_base=_normalize_grid(base)
        aug=augment_with_corridor_waves(norm_base)
        ratio=_ff_count(aug)/max(1,_omega_open(aug))
        score=(ratio, -len(aug), _ff_count(aug))
        if best is None or score>best[0]:
            best=(score,aug)
    return best[1]

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()