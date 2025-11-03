# EVOLVE-BLOCK-START

from math import gcd

# --- Core utilities ---

def normalize_grid(intervals):
    if not intervals:
        return []
    pts = sorted({x for iv in intervals for x in iv})
    mp = {v: i*2 for i, v in enumerate(pts)}
    norm = [(mp[l], mp[r]) for (l, r) in intervals]
    g = 0
    for l, r in norm:
        g = gcd(g, abs(l)); g = gcd(g, abs(r))
    if g > 1:
        norm = [(l//g, r//g) for l, r in norm]
    return norm

def firstfit(intervals):
    last = []
    for l, r in intervals:
        placed = False
        for i in range(len(last)):
            if l >= last[i]:
                last[i] = r
                placed = True
                break
        if not placed:
            last.append(r)
    return len(last)

def omega(intervals):
    ev = []
    for l, r in intervals:
        if l<r:
            ev.append((l,1)); ev.append((r,-1))
    ev.sort(key=lambda e:(e[0], 0 if e[1]==-1 else 1))
    cur = best = 0
    for _, t in ev:
        cur += t
        if cur>best: best=cur
    return best

# --- Pattern builder with tiling‐cycle and interleaving ---

OFFSETS_A = (2,6,10,14)
OFFSETS_B = (1,5,9,13)
OFFSETS_C = (3,7,11,15)
OFFSETS_D = (0,4,8,12)
PATTERNS   = [OFFSETS_A, OFFSETS_B, OFFSETS_C, OFFSETS_D]

BLOCKERS    = [(1,5),(12,16),(4,9),(8,13)]

def build_tiling_cycle(depth, extra_first=False, extra_last=False, schedule='interleaved'):
    T = [(0.0,1.0)]
    for lvl in range(depth):
        lo = min(l for l,r in T)
        hi = max(r for l,r in T)
        delta = hi-lo
        pat = PATTERNS[lvl % 4]
        offs = list(pat)
        if extra_first and lvl==0:
            offs.append(max(pat)+4)
        if extra_last and lvl==depth-1:
            offs.append(max(pat)+4)
        # precompute scaled blockers
        blks = [(delta*a, delta*b) for (a,b) in BLOCKERS]
        # interleaved insertion
        S = []
        if schedule=='interleaved':
            for i, start in enumerate(offs):
                off = delta*start - lo
                for (l,r) in T:
                    S.append((l+off, r+off))
                # cycle through blockers too
                S.append(blks[i % len(blks)])
        else:
            # split or after
            if schedule=='split':
                h = len(offs)//2
                for start in offs[:h]:
                    off = delta*start - lo
                    for (l,r) in T: S.append((l+off,r+off))
                S.extend(blks)
                for start in offs[h:]:
                    off = delta*start - lo
                    for (l,r) in T: S.append((l+off,r+off))
            else:  # after
                for start in offs:
                    off = delta*start - lo
                    for (l,r) in T: S.append((l+off,r+off))
                S.extend(blks)
        T = S
    return T

# --- Pruning ---

def prune(intervals, target_ratio):
    cur = list(intervals)
    norm = normalize_grid(cur)
    if not norm: return cur
    base_cols = firstfit(norm)
    base_om   = omega(norm) or 1
    if target_ratio is None:
        target_ratio = base_cols/base_om
    def length(iv): return iv[1]-iv[0]
    changed = True
    while changed:
        changed = False
        order = sorted(range(len(cur)), key=lambda i:(-length(cur[i]),i))
        for i in order:
            cand = cur[:i]+cur[i+1:]
            cn = normalize_grid(cand)
            if not cn: continue
            ccols = firstfit(cn)
            com   = omega(cn) or 1
            if ccols/com >= target_ratio-1e-12:
                cur = cand; changed=True; break
    return cur

# --- Search and selection ---

def construct_intervals():
    best = None  # (ratio, -n, raw)
    # try depths 4 and 5
    for depth in (4,5):
        for extra in (False, True):
            for extra2 in (False, True):
                for sched in ('interleaved','split','after'):
                    T = build_tiling_cycle(depth, extra, extra2, sched)
                    N = normalize_grid(T)
                    if not N: continue
                    om = omega(N); 
                    if om<=0: continue
                    cols = firstfit(N)
                    ratio = cols/om
                    key = (ratio, -len(T))
                    if best is None or key>best[0]:
                        best = (key, T, ratio)
    if best is None:
        # fallback simple
        T = build_tiling_cycle(4, False, False, 'after')
        return normalize_grid(T)
    raw = best[1]
    target = best[2]
    pr = prune(raw, target)
    return normalize_grid(pr)

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()