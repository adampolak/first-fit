# EVOLVE-BLOCK-START

from itertools import product

def overlaps(a, b):
    # open‐interval overlap
    (l1,r1),(l2,r2)=a,b
    return max(l1,l2)<min(r1,r2)

def firstfit_colors(intervals):
    # Greedy FirstFit via tracking last end per color class
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
    # Sweep‐line to compute max overlap (ω) for open intervals
    events=[]
    for l,r in intervals:
        if l<r:
            events.append((l,1))
            events.append((r,-1))
    events.sort(key=lambda e:(e[0],0 if e[1]==-1 else 1))
    cur=best=0
    for _,t in events:
        cur+=t
        if cur>best:
            best=cur
    return best

def normalize_to_grid(intervals):
    # Map unique endpoints to even integers preserving order
    pts=sorted(set(x for l,r in intervals for x in (l,r)))
    coord={e:2*i for i,e in enumerate(pts)}
    return [(coord[l],coord[r]) for l,r in intervals]

def build_pattern(p, depth, start0=2, step=4, connectors=None):
    """
    Build a p‐copy recursive adversary:
      - p: branching factor
      - depth: recursion levels
      - connectors: list of (a,b) multiplier pairs for long intervals
    """
    T=[(0.0,1.0)]
    for _ in range(depth):
        lo=min(l for l,r in T)
        hi=max(r for l,r in T)
        delta=hi-lo
        # generate p translated/scaled copies
        offsets=[start0+step*i for i in range(p)]
        S=[]
        for m in offsets:
            off=delta*m - lo
            for l,r in T:
                S.append((l+off, r+off))
        # append connector intervals
        if connectors:
            for a,b in connectors:
                S.append((delta*a - lo, delta*b - lo))
        T=S
    return T

def construct_intervals():
    """
    Enumerate patterns (p, depth, connector template), choose the one
    maximizing FirstFit/ω under a size bound, then normalize.
    """
    best=None
    max_n=2000
    # try p∈{3,4,5}
    for p in (3,4,5):
        # define two connector styles
        offsets=[2+4*i for i in range(p)]
        conn_adj=[(offsets[i]-1, offsets[i+1]+1) for i in range(p-1)]
        conn_glob=[(offsets[0]-1, offsets[-1]+1)]
        for connectors in (conn_adj, conn_glob, conn_adj+conn_glob):
            for depth in (3,4,5):
                T=build_pattern(p, depth, connectors=connectors)
                Tn=normalize_to_grid(T)
                n=len(Tn)
                if n==0 or n>max_n:
                    continue
                om=clique_number(Tn)
                if om<=0:
                    continue
                cols=firstfit_colors(Tn)
                ratio=cols/om
                # tuple ordering: maximize ratio, then minimize size, then maximize colors
                cand=(ratio, -n, cols, om, T)
                if best is None or cand>best:
                    best=cand
    # fallback to trivial single‐interval if search failed
    if best is None:
        T=[(0.0,1.0)]
    else:
        T=best[4]
    return normalize_to_grid(T)

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()