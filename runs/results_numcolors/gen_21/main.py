# EVOLVE-BLOCK-START

def construct_intervals(depth=4):
    """
    Recursive construction of interval sequence inspired by Kierstead-Trotter style
    blow-up (see Figure 4 in https://arxiv.org/abs/1506.00192), enhanced with:
      - alternating clone order (reverse every other block)
      - short preload blockers inside each cloned block
      - short bridges between adjacent blocks
    depth: number of recursive expansions (default 4)
    Returns list of (l,r) intervals in presentation order.
    """
    # Base pattern T_0
    T = [(0.0, 1.0)]
    # Parameters to control added gadgets (kept small to avoid large omega)
    preload_blockers = 2      # number of short blockers prepended in each block
    blocker_length = 0.18     # length of each preload blocker
    blocker_gap = 0.06        # spacing between preload blockers inside a block
    for level in range(depth):
        # compute current span
        lo = min(l for l, r in T)
        hi = max(r for l, r in T)
        span = hi - lo
        # make 4 scaled copies at offsets 2,6,10,14 times span (core pattern)
        S = []
        offsets = [2, 6, 10, 14]
        # alternate clone order every other level to disrupt FF reuse
        reverse_flag = (level % 2 == 1)
        for idx, off in enumerate(offsets):
            base = off * span - lo
            # insert small preload blockers inside the block before clones
            for b in range(preload_blockers):
                start = base + (0.08 + b * (blocker_length + blocker_gap)) * span / span
                # blockers are tiny intervals placed inside the block area [off*span, (off+1)*span]
                S.append((off * span + 0.08 * span + b * (blocker_length + blocker_gap) * span,
                          off * span + 0.08 * span + b * (blocker_length + blocker_gap) * span + blocker_length * span))
            # append cloned copy of T (possibly reversed)
            if reverse_flag:
                for (l, r) in reversed(T):
                    S.append((l + base, r + base))
            else:
                for (l, r) in T:
                    S.append((l + base, r + base))
            # short bridge to next block (couples color choices lightly)
            if idx + 1 < len(offsets):
                next_off = offsets[idx + 1]
                # bridge spans a small portion that overlaps the tail of this block and the head of next
                S.append(( (off + 0.55) * span, (next_off + 0.35) * span ))
        # Add the 4 connector intervals as in Figure 4 (keeps global coupling)
        S.append((1 * span, 5 * span))    # connects first two blocks
        S.append((12 * span, 16 * span))  # connects last two blocks
        S.append((4 * span, 9 * span))    # cross connector
        S.append((8 * span, 13 * span))   # cross connector
        # keep T as floats for more varied fractional overlaps
        T = S
    return T

def run_experiment(**kwargs):
    """Main called by evaluator"""
    return construct_intervals()

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()