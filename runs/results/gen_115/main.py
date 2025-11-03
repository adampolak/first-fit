# EVOLVE-BLOCK-START

def construct_intervals():
    """
    Build a 5‐branch fractal with rich connectors, zip interleaving, and alternating reversal.
    Returns a normalized list of integer intervals.
    """
    # Base seed: one unit interval
    T = [(0.0, 1.0)]
    # Parameters
    depth = 4
    offsets = [2, 6, 10, 14, 18]      # 5 copies per level
    # Connect adjacent copies plus an end‐cap connector
    connectors = [(offsets[i]-1, offsets[i+1]+1) for i in range(len(offsets)-1)]
    connectors.append((offsets[0]-2, offsets[-1]+2))

    for level in range(depth):
        lo = min(l for l, r in T)
        hi = max(r for l, r in T)
        delta = hi - lo
        center = (lo + hi) / 2.0

        # Build each copy, alternating reversed order
        copy_lists = []
        for idx, start in enumerate(offsets):
            off = delta * start - center
            seq = T if (idx % 2 == 0) else list(reversed(T))
            copy_lists.append([(l + off, r + off) for l, r in seq])

        # Zip‐interleave all copies
        S = []
        block_size = len(copy_lists[0])
        for i in range(block_size):
            for lst in copy_lists:
                S.append(lst[i])

        # Add all blockers (connectors)
        for a, b in connectors:
            S.append((delta * a - center, delta * b - center))

        # Next iteration uses S
        T = S

    # Normalize to small integer grid (even spacing)
    pts = sorted(set(x for seg in T for x in seg))
    coord = {v: 2*i for i, v in enumerate(pts)}
    normalized = [(coord[l], coord[r]) for (l, r) in T]
    return normalized

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()