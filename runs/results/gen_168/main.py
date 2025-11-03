# EVOLVE-BLOCK-START

from functools import lru_cache

def _normalize_grid(intervals):
    endpoints = sorted({x for seg in intervals for x in seg})
    coord = {v:2*i for i,v in enumerate(endpoints)}
    return [(coord[l],coord[r]) for (l,r) in intervals]

def _ff_count(intervals):
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

def _omega_open(intervals):
    ev=[]
    for l,r in intervals:
        if l<r:
            ev.append((l,1)); ev.append((r,-1))
    ev.sort(key=lambda x:(x[0],0 if x[1]==-1 else 1))
    cur=best=0
    for _,t in ev:
        cur+=t
        best=max(best,cur)
    return best

# four offset patterns A,B,C,D
PATTERNS=[(2,6,10,14),(1,5,9,13),(3,7,11,15),(0,4,8,12)]
# canonical blockers
BLOCKERS=[(1,5),(12,16),(4,9),(8,13)]

@lru_cache(maxsize=None)
def _build(level, pattern_idx):
    """
    Recursively build at given level using pattern PATTERNS[pattern_idx %4].
    Inserts waves: short if level%2==0, long if level%3==0.
    """
    if level==0:
        return ((0.0,1.0),)
    prev=_build(level-1,pattern_idx-1)
    lo=min(l for l,r in prev); hi=max(r for l,r in prev)
    delta=hi-lo
    center=(lo+hi)/2.0
    offs=PATTERNS[pattern_idx%4]
    # place copies
    S=[]
    for st in offs:
        off=delta*st - lo
        for l,r in prev:
            S.append((l+off,r+off))
    # blockers
    for a,b in BLOCKERS:
        S.append((delta*a,delta*b))
    # short wave
    if level%2==0:
        S.append((delta*3,delta*5))
        S.append((delta*5,delta*7))
    # long wave
    if level%3==0:
        S.append((delta*7,delta*11))
        S.append((delta*11,delta*17))
    return tuple(S)

def _prune(intervals):
    """
    Single-pass left-to-right removal if ratio unchanged.
    """
    cur=list(intervals)
    # normalize for evaluation
    def eval_ratio(T):
        Tn=_normalize_grid(T)
        om=_omega_open(Tn)
        if om==0: return -1
        return _ff_count(Tn)/om
    base=eval_ratio(cur)
    i=0
    while i<len(cur):
        cand=cur[:i]+cur[i+1:]
        if eval_ratio(cand)>=base:
            cur=cand
            # do not increment i to re-examine at same index
        else:
            i+=1
    return cur

def construct_intervals(iterations=4, normalize=True):
    """
    Build at fixed depth=4 using cyclic patterns and waves, then prune.
    """
    DEPTH=4
    T= list(_build(DEPTH,DEPTH))
    T_pruned=_prune(T)
    return _normalize_grid(T_pruned) if normalize else T_pruned

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()