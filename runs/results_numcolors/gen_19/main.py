# EVOLVE-BLOCK-START

def construct_intervals():
    """
    Two-phase deterministic construction with preload blockers, multiple seeds,
    and occasional long sentries to raise FirstFit pressure while keeping omega small.

    Strategy highlights:
      - Seed with multiple disjoint unit intervals to increase early contention.
      - Two-phase iteration: first phase uses the classical 4-block offsets with
        preload blockers; second phase switches to a different start template
        to create new overlap patterns.
      - Small preload blockers inside each block force low-color occupation early.
      - A handful of long sentry intervals (added once at the beginning) reserve
        low colors for a long span and push FirstFit upward.
      - Alternating clone orders, rotated block order, and short bridges between
        neighboring blocks further disturb FirstFit reuse.
    """
    # overall recursion depth (kept similar to prior; total size remains controlled)
    depth = 6

    # Base seed: multiple disjoint unit intervals (increases early overlaps without
    # forming large cliques because they are disjoint)
    T = [(0.0, 1.0), (2.0, 3.0), (4.0, 5.0)]

    # Start-pattern cycles used in phase 1 (classical offsets) and phase 2 (shifted)
    start_patterns_phase1 = [
        (2, 6, 10, 14),
        (1, 5, 9, 13),
        (3, 7, 11, 15),
        (2, 4, 8, 12),
    ]
    start_patterns_phase2 = [
        (2, 6, 10, 14, 18),   # slightly more blocks but still spaced
        (1, 5, 9, 13, 17),
        (3, 7, 11, 15, 19),
        (0, 4, 8, 12, 16),
    ]

    # Two small template banks for the two phases (keeps connector count modest)
    template_bank_p1 = [
        ((1, 5), (12, 16), (4, 9), (8, 13)),
        ((0, 4), (5, 9), (2, 6), (7, 11)),
    ]
    template_bank_p2 = [
        ((1, 6), (10, 15), (3, 8), (12, 17)),
        ((0, 5), (7, 12), (2, 7), (9, 14)),
    ]

    # Preload blocker parameters (short intervals placed near the start of each block)
    preload_count = 2
    preload_len = 0.18  # length of each small blocker (in block-span units)
    preload_gap = 0.06  # small step between preload blockers inside a block

    # Add a small set of long sentries once to occupy low colors early
    # (these are long spans relative to T's current span but they are added only once)
    added_sentries = False

    for it in range(depth):
        lo = min(l for l, r in T)
        hi = max(r for l, r in T)
        span = hi - lo

        # Decide which phase we're in; keep first 4 iterations as phase1, rest as phase2
        if it < 4:
            starts = list(start_patterns_phase1[it % len(start_patterns_phase1)])
            template_bank = template_bank_p1
        else:
            starts = list(start_patterns_phase2[it % len(start_patterns_phase2)])
            template_bank = template_bank_p2

        # vary block order: rotate / reverse to change arrival interleaving
        if it % 4 == 0:
            block_order = starts
        elif it % 4 == 1:
            block_order = list(reversed(starts))
        elif it % 4 == 2:
            block_order = starts[1:] + starts[:1]
        else:
            # alternate every other time with a 2-shift
            block_order = starts[2:] + starts[:2]

        S = []

        # Add sentries only once at the very beginning to occupy early small colors
        if not added_sentries:
            # two long sentries with modest overlap (avoid raising omega too much)
            S.append((span * 0.5, span * 2.5))
            S.append((span * 1.2, span * 3.2))
            added_sentries = True

        for idx, st in enumerate(block_order):
            # Preload blockers: small intervals placed early within the block to
            # occupy low-numbered colors locally.
            base_block_left = span * st
            for b in range(preload_count):
                left = base_block_left + preload_gap * b * span
                right = left + preload_len * span
                S.append((left, right))

            # choose clone order to alternate symmetry
            if (it + idx) % 2 == 0:
                clone_order = T
            else:
                clone_order = list(reversed(T))

            # append cloned copies offset by block start
            base_off = span * st - lo
            for (l, r) in clone_order:
                S.append((l + base_off, r + base_off))

            # short bridge to the next block to couple colors (but short enough
            # to avoid creating large cliques)
            if idx + 1 < len(block_order):
                nxt = block_order[idx + 1]
                left = span * (st + 0.55)
                right = span * (nxt + 0.25)
                S.append((left, right))

        # append a modest set of connectors from the chosen template bank
        template = template_bank[it % len(template_bank)]
        for (a, b) in template:
            S.append((span * a, span * b))

        # occasionally add a cross connector to mix older layers and new ones
        if it % 2 == 0:
            S.append((span * 4.5, span * 9.5))

        # update T for next iteration
        T = S

    return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()