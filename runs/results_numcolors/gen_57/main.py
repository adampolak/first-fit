# EVOLVE-BLOCK-START

def construct_intervals(rounds=3, seed_lo=0.0, seed_scale=1.0):
    """
    Enhanced deterministic rotating multiscale construction with interleaved blocks.

    Improvements vs. prior version:
    - Default rounds increased to 3 for deeper recursion.
    - Build each round as separate blocks (one per start) and present them
      in a round-robin interleaved order to force many active colors.
    - Add thin bridges after interleaving so they overlap many active items.
    - Vary tower span and layers per round to create staggered stacks.
    - Slightly larger micro-gadget and a couple more caps to increase local pressure.
    """

    # -- Key knobs (exposed for experimentation) --
    phase2_rounds = 1           # small fine phase to couple scales
    bridge_fraction = 0.38      # fraction controlling internal bridge placement
    caps_per_round = 3          # a few more caps for stronger local pressure
    micro_scale = 0.22          # micro-gadget scale inside a mid block
    interleave_blocks = True    # present blocks in interleaved (round-robin) order
    # ------------------------------------------------

    # Seed: 4 disjoint unit intervals to provide richer early interactions
    T = [
        (seed_lo + 0.0 * seed_scale, seed_lo + 1.0 * seed_scale),
        (seed_lo + 2.0 * seed_scale, seed_lo + 3.0 * seed_scale),
        (seed_lo + 4.0 * seed_scale, seed_lo + 5.0 * seed_scale),
        (seed_lo + 6.0 * seed_scale, seed_lo + 7.0 * seed_scale),
    ]

    # Rotating cycle of start patterns (now include 3- and 4-start variants)
    start_patterns = [
        (2, 6, 10, 14),
        (1, 5, 9, 13),
        (3, 7, 11, 15),
        (2, 4, 8, 12),
    ]

    # Template bank: 4-interval gadget variants (keeps omega small)
    template_bank = [
        ((1.0, 5.0), (12.0, 16.0), (4.0, 9.0), (8.0, 13.0)),
        ((0.5, 4.5), (11.0, 15.0), (3.5, 8.5), (7.0, 12.0)),
        ((1.0, 4.0), (6.0, 9.0), (3.0, 7.0), (9.0, 13.0)),
        ((2.0, 6.0), (7.0, 11.0), (0.0, 3.0), (10.0, 14.0)),
    ]

    for it in range(max(1, int(rounds))):
        lo = min(l for l, r in T)
        hi = max(r for l, r in T)
        delta = hi - lo if hi > lo else 1.0

        # choose rotating start pattern and possibly truncate to 3 or 4 starts
        starts = list(start_patterns[it % len(start_patterns)])
        # sometimes use first 3 starts to vary block counts deterministically
        if it % 2 == 0:
            starts = starts[:3]

        # permute block order deterministically to disrupt reuse patterns
        if it % 4 == 0:
            block_order = starts
        elif it % 4 == 1:
            block_order = list(reversed(starts))
        elif it % 4 == 2:
            block_order = starts[1:] + starts[:1]
        else:
            # rotate-right
            block_order = starts[-1:] + starts[:-1]

        # -- Build blocks (store each translated clone as its own list) --
        blocks = []
        for bi, st in enumerate(block_order):
            clone_order = T if ((it + bi) % 2 == 0) else list(reversed(T))
            shift = delta * st - lo
            block = [(l + shift, r + shift) for (l, r) in clone_order]
            blocks.append(block)

        # -- Interleave blocks (round-robin) to maximize simultaneous active colors --
        S = []
        if interleave_blocks and len(blocks) > 1:
            maxlen = max(len(b) for b in blocks)
            order = list(range(len(blocks)))
            # occasionally reverse the visitation order to diversify overlaps
            if it % 2 == 1:
                order = order[::-1]
            for i in range(maxlen):
                for idx in order:
                    blk = blocks[idx]
                    if i < len(blk):
                        S.append(blk[i])
        else:
            # simple concatenation (fallback)
            for blk in blocks:
                S.extend(blk)

        # -- Add thin bridges AFTER interleaving so they meet many active colors --
        for idx in range(len(block_order) - 1):
            st = block_order[idx]
            nxt = block_order[idx + 1]
            # place bridge interior-to-interior with short length
            b_lo = delta * (st + bridge_fraction)
            # small bridge that ends inside next block's interior region
            b_len = max(0.08 * delta, 0.4)
            S.append((b_lo, b_lo + b_len))

        # -- add rotating template gadget at global scale --
        template = template_bank[it % len(template_bank)]
        for (a, b) in template:
            S.append((delta * a, delta * b))

        # -- micro gadget inside the middle block (small scale) --
        mid_start = block_order[len(block_order) // 2]
        for (a, b) in template:
            S.append((delta * (mid_start + micro_scale * a),
                      delta * (mid_start + micro_scale * b)))

        # -- Towers: vary layers and span per round to create staggered coupling --
        tower_layers = 2 + (it % 3)  # 2..4 layers across rounds
        # alternate span to avoid creating single-point cliques
        span_blocks = 1.6 if (it % 2 == 0) else 2.05
        for bi, st in enumerate(block_order[:-1]):
            for layer in range(tower_layers):
                off = layer * (0.42 / max(1, tower_layers - 1))
                t_lo = delta * (st + 0.18 + off)
                t_hi = delta * (st + span_blocks + 0.18 + off)
                if t_hi > t_lo:
                    S.append((t_lo, t_hi))

        # -- small caps to nudge FirstFit (spread positions deterministically) --
        max_start = max(starts) if starts else 1
        cap_positions = []
        for i in range(caps_per_round):
            # choose fractional positions to hit different relative phases
            pos = 1 + int((max_start - 1) * (i + 1) / (caps_per_round + 1))
            cap_positions.append(pos)
        for j in cap_positions:
            cap_lo = delta * (j + 0.12)
            cap_hi = cap_lo + max(0.06 * delta, 0.35)
            S.append((cap_lo, cap_hi))

        # prepare for next round
        T = S

    # -- Fine-scale second phase to couple scales --
    lo = min(l for l, r in T)
    hi = max(r for l, r in T)
    delta = hi - lo if hi > lo else 1.0
    delta2 = max(1.0, delta / 4.0)  # slightly smaller fine scale for denser interactions

    fine_starts_cycle = [
        (1, 4, 7),
    ]

    for p in range(phase2_rounds):
        starts = list(fine_starts_cycle[p % len(fine_starts_cycle)])
        blocks = []
        for bi, st in enumerate(starts):
            clone_order = T if ((p + bi) % 2 == 0) else list(reversed(T))
            shift = delta2 * st - lo
            block = [(l * (delta2 / delta) + shift, r * (delta2 / delta) + shift) for (l, r) in clone_order]
            blocks.append(block)

        # interleave fine blocks too
        S = []
        maxlen = max(len(b) for b in blocks)
        for i in range(maxlen):
            for blk in blocks:
                if i < len(blk):
                    S.append(blk[i])

        # fine bridges
        for idx in range(len(starts) - 1):
            st = starts[idx]
            b_lo = delta2 * (st + 0.28)
            b_hi = b_lo + max(0.03 * delta2, 0.15)
            S.append((b_lo, b_hi))

        # a small cap-coupler
        S.append((delta2 * (starts[0] + 0.35), delta2 * (starts[-1] - 0.35)))

        T = S

    return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()