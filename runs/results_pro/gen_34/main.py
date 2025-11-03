# EVOLVE-BLOCK-START

def construct_intervals():
    import math

    # --- Initial seeds -----------------------------------------------------
    T = [(0.0, 1.0), (3.0, 4.0)]

    # --- Pattern banks -----------------------------------------------------
    START_PATTERNS = [
        (2, 6, 10, 14),
        (3, 7, 11, 15),
        (4, 8, 12, 16),
        (5, 9, 13, 17),
        (6, 10, 14, 18),
    ]
    BRIDGE_SETS = [
        [(1, 5),  (4, 9),  (8, 13), (12, 16)],
        [(2, 6),  (5, 10), (9, 14), (13, 17)],
        [(3, 7),  (6, 11), (10, 15),(14, 18)],
        [(4, 8),  (7, 12), (11, 16),(15, 19)],
        [(5, 9),  (8, 13), (12, 17),(16, 20)],
    ]
    GAMMA = [1.00, 1.20, 0.80, 1.10, 1.05]
    JITTER = [0.0, 0.02, -0.015, 0.04, -0.03]

    # --- Dynamic depth schedule & recursive expansion ----------------------
    DEPTH_SEQ = [3, 4, 5, 4]
    iter_idx = 0
    for depth in DEPTH_SEQ:
        for _ in range(depth):
            lo = min(l for l, r in T)
            hi = max(r for l, r in T)
            span = hi - lo if hi > lo else 1.0

            gamma = GAMMA[iter_idx % len(GAMMA)]
            delta = span * gamma
            jit = JITTER[iter_idx % len(JITTER)] * span * 0.005

            pattern = START_PATTERNS[iter_idx % len(START_PATTERNS)]
            bridges = BRIDGE_SETS[iter_idx % len(BRIDGE_SETS)]

            S = []
            # replicate with jittered offsets
            for s in pattern:
                offset = delta * s - lo + jit
                for (l, r) in T:
                    S.append((l + offset, r + offset))
            # add jittered bridge intervals
            for (a, b) in bridges:
                S.append((delta * a + jit, delta * b + jit))

            T = S
            iter_idx += 1

    # --- FirstFit utilities ------------------------------------------------
    def overlap(a, b):
        return a[0] < b[1] and b[0] < a[1]

    def firstfit(ints):
        cols = []
        for i, iv in enumerate(ints):
            used = {cols[j] for j in range(i) if overlap(iv, ints[j])}
            c = 1
            while c in used:
                c += 1
            cols.append(c)
        return cols

    # --- Targeted capping phase ---------------------------------------------
    cols = firstfit(T)
    maxcol = max(cols) if cols else 0

    # Add one short "cap" interval per existing color to force new FF colors
    lo_all = min(l for l, r in T)
    hi_all = max(r for l, r in T)
    total_span = hi_all - lo_all if hi_all > lo_all else 1.0

    caps = []
    for color in range(1, min(maxcol, 12) + 1):
        centers = [
            (T[i][0] + T[i][1]) / 2.0
            for i, c in enumerate(cols) if c == color
        ]
        if not centers:
            continue
        left, right = min(centers), max(centers)
        length = max((right - left) * 0.2, total_span * 0.01)
        mid = (left + right) / 2.0
        pad = total_span * 0.005
        cap = (mid - length / 2 + pad, mid + length / 2 + pad)
        caps.append(cap)

    T += caps

    # --- Normalize to integer grid -----------------------------------------
    lo = min(l for l, r in T)
    if lo < 0:
        shift = -lo + 1.0
        T = [(l + shift, r + shift) for (l, r) in T]

    hi = max(r for l, r in T)
    scale = max(1.0, 1000.0 / hi)

    out = []
    for (l, r) in T:
        L = int(math.floor(l * scale))
        R = int(math.ceil(r * scale))
        if R <= L:
            R = L + 1
        out.append((L, R))

    # deduplicate exactly-equal intervals
    seen = set()
    final = []
    for iv in out:
        if iv not in seen:
            seen.add(iv)
            final.append(iv)

    return final

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()