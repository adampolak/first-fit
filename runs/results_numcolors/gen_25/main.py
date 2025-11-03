# EVOLVE-BLOCK-START

def construct_intervals():
    """
    Hybrid recursive interval generator combining Kierstead–Trotter blow-up
    with arrival-order heuristics for FirstFit pressure.
    Returns list of integer (l, r) intervals.
    """
    depth = 4             # recursion depth
    branching = 4         # number of scaled copies per level
    preload = 2           # number of small blockers per block
    caps_per_round = 3    # number of sparse caps each round

    # base pattern
    T = [(0, 1)]
    for round_idx in range(depth):
        lo = min(l for l, r in T)
        hi = max(r for l, r in T)
        span = hi - lo

        S = []
        # block start offsets as in Figure 4
        starts = [2 + 4 * i for i in range(branching)]

        for idx, start in enumerate(starts):
            base = start * span - lo

            # preload blockers: small intervals to occupy low colors
            small_len = span * 0.1
            for j in range(preload):
                l0 = (start + 0.1 + j * 0.15) * span - lo
                S.append((int(l0), int(l0 + small_len)))

            # clone or reversed‐clone of previous round
            if idx % 2 == 1:
                for (l, r) in reversed(T):
                    S.append((int(l + base), int(r + base)))
            else:
                for (l, r) in T:
                    S.append((int(l + base), int(r + base)))

            # short bridge to next block
            if idx + 1 < len(starts):
                next_start = starts[idx + 1]
                l1 = (start + 0.5) * span - lo
                r1 = (next_start + 0.3) * span - lo
                S.append((int(l1), int(r1)))

        # adjacent connectors (propagate FF usage between neighbors)
        for i in range(branching - 1):
            a = (1 + 4 * i) * span
            b = (5 + 4 * i) * span
            S.append((int(a), int(b)))

        # the four global connectors from Figure 4
        S.append((int(1 * span),  int(5 * span)))
        S.append((int(12 * span), int(16 * span)))
        S.append((int(4 * span),  int(9 * span)))
        S.append((int(8 * span),  int(13 * span)))

        # sparse caps to raise FF colors without enlarging clique
        for c in range(1, caps_per_round + 1):
            j = 6 * c
            l2 = (j - 1.5) * span
            r2 = (j + 1.5) * span
            S.append((int(l2), int(r2)))

        # two long sentries in the first round only
        if round_idx == 0:
            S.append((0,                int(3 * span)))
            S.append((int(10 * span),   int(14 * span)))

        # prepare for next round
        T = S

    return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()