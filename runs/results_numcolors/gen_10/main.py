# EVOLVE-BLOCK-START

def construct_intervals(depth=4):
    """
    Recursive construction of interval sequence inspired by Kierstead-Trotter style
    blow-up (see Figure 4 in https://arxiv.org/abs/1506.00192), with one extra
    level of recursion to push the FirstFit ratio beyond the previous 2.33.
    depth: number of recursive expansions (default 3)
    Returns list of (l,r) intervals in presentation order.
    """
    # Base pattern T_0
    T = [(0.0, 1.0)]
    for _ in range(depth):
        # compute current span
        lo = min(l for l, r in T)
        hi = max(r for l, r in T)
        span = hi - lo
        # we will make 4 scaled copies at offsets 2,6,10,14 times span
        S = []
        offsets = [2, 6, 10, 14]
        for idx, off in enumerate(offsets):
            base = off * span - lo
            # replicate block
            for (l, r) in T:
                S.append((l + base, r + base))
            # ladder connector to next block
            if idx < len(offsets) - 1:
                next_off = offsets[idx + 1]
                S.append(((off + 1) * span, (next_off - 1) * span))
        # cross-layer connectors to couple nonconsecutive blocks
        # connects block0 to block2
        S.append(((offsets[0] + 3) * span, (offsets[2] + 7) * span))
        # connects block1 to block3
        S.append(((offsets[1] + 3) * span, (offsets[3] + 7) * span))
        # original figure connectors to maintain structure
        S.append((1 * span, 5 * span))
        S.append((12 * span, 16 * span))
        T = S
    return T

def run_experiment(**kwargs):
    return construct_intervals(depth=4)

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()