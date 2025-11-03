# EVOLVE-BLOCK-START

from math import gcd

# --- overlap, FirstFit, clique on floats ---

def overlaps(a,b):
    (l1,r1),(l2,r2)=a,b
    return max(l1,l2)<min(r1,r2)

def firstfit_colors(T):
    colors=[]
    for iv in T:
        placed=False
        for c in colors:
            if not any(overlaps(iv,u) for u in c):
                c.append(iv)
                placed=True
                break
        if not placed:
            colors.append([iv])
    return len(colors)

def clique_number(T):
    ev=[]
    for l,r in T:
        if l<r:
            ev.append((l,1))
            ev.append((r,-1))
    ev.sort(key=lambda x:(x[0], x[1]))
    cur=best=0
    for _,t in ev:
        cur+=t
        if cur>best: best=cur
    return best

# --- strict pruning on floats (slightly relaxed: allow removals that preserve or strengthen FF) ---

def prune_strict_floats(intervals, target_cols, target_om):
    cur=list(intervals)
    # sort by descending length to try to remove long redundant ones first
    def length(iv): return iv[1]-iv[0]
    changed=True
    while changed:
        changed=False
        order=sorted(range(len(cur)), key=lambda i:(-length(cur[i]),i))
        for i in order:
            cand=cur[:i]+cur[i+1:]
            # accept removal only if clique unchanged and FF not decreased below target
            if clique_number(cand)==target_om and firstfit_colors(cand)>=target_cols:
                cur=cand
                changed=True
                break
    return cur

# --- normalize floats -> small ints ---

def normalize_intervals(T):
    if not T: return []
    pts=sorted({x for seg in T for x in seg})
    # map to even integers
    mp={}
    c=0
    for x in pts:
        mp[x]=c; c+=2
    L=[(mp[l],mp[r]) for l,r in T]
    # shift min to 0
    mn=min(min(a,b) for a,b in L)
    L=[(a-mn,b-mn) for a,b in L]
    # divide out gcd
    g=0
    for a,b in L:
        g=gcd(g,abs(a)); g=gcd(g,abs(b))
    if g>1:
        L=[(a//g,b//g) for a,b in L]
    return L

# --- interleaving helpers ---

def zip_interleave(list_of_lists):
    if not list_of_lists:
        return []
    m = max(len(lst) for lst in list_of_lists)
    S=[]
    for i in range(m):
        for lst in list_of_lists:
            if i < len(lst):
                S.append(lst[i])
    return S

# --- build one round of 4-copy + 4-blockers with adversarial ordering ---

def expand_once(T, level):
    lo=min(l for l,r in T)
    hi=max(r for l,r in T)
    delta=hi-lo
    offsets=(2,6,10,14)
    # build copies as separate lists (so we can interleave)
    copies=[]
    for idx,s in enumerate(offsets):
        off=delta*s - lo
        seq=[(l+off,r+off) for l,r in T]
        # occasionally reverse every other copy to break alignment
        if (level % 2)==1 and (idx%2)==1:
            seq=list(reversed(seq))
        copies.append(seq)
    # choose interleaving style by level: alternate 'zip' and 'block'
    if (level % 2)==1:
        S_copies = zip_interleave(copies)
    else:
        S_copies = []
        for lst in copies:
            S_copies.extend(lst)
    # blockers remain long intervals appended at end (helps keep omega small)
    blockers=[(1,5),(12,16),(4,9),(8,13)]
    S = S_copies + [(delta*a, delta*b) for a,b in blockers]
    return S

# --- main construction with per-level pruning and improved base seed ---

def construct_intervals():
    # start with a slightly richer base seed: two overlapping short intervals
    T=[(0.0,1.0),(0.5,1.5)]
    # depth 4 is a sweet spot in experiments; keep but use adversarial interleaving
    for i in range(4):
        S=expand_once(T, i)
        # compute local targets on floats
        target_om=clique_number(S)
        target_cols=firstfit_colors(S)
        # prune strictly on floats but allow removals that preserve or strengthen FF
        T=prune_strict_floats(S, target_cols, target_om)
    # normalize to small ints
    return normalize_intervals(T)

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()