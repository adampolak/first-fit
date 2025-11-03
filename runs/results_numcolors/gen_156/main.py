# EVOLVE-BLOCK-START

def construct_intervals(rounds=6,
                        rotate_starts=True,
                        reverse_block_parity=True,
                        interleave_blocks=True,
                        phase2_iters=2):
    """
    Improved KT-style spine with parity interleaving, spine multiplier, and cross4 connectors.

    Signature maintained for compatibility. Returns a list of (l, r) integer tuples
    representing open intervals in the order presented to FirstFit.
    """

    # Hard capacity guard to keep total count comfortably < 10k
    CAP = 9800

    # Four core strong templates (T1-T4). We rotate among them round-robin.
    template_bank = [
        [2, 6, 10, 14],  # classic KT (T1)
        [1, 5, 9, 13],   # left-shifted (T2)
        [3, 7, 11, 15],  # right-shifted (T3)
        [4, 8, 12, 16],  # stretched-right (T4)
    ]

    # Single unit seed (keeps omega low early)
    T = [(0, 1)]

    # Estimate how many full spine rounds we can do without exceeding CAP
    def max_rounds_within_cap(initial_size, max_rounds):
        sz = initial_size
        done = 0
        for _ in range(max(0, int(max_rounds))):
            nxt = 4 * sz + 4
            if nxt > CAP:
                break
            sz = nxt
            done += 1
        return done, sz

    depth, _ = max_rounds_within_cap(len(T), rounds)

    # Helper to build translated blocks (with optional internal reversal)
    def build_blocks(current_T, starts, delta, K=1, reverse_block_parity=False, round_idx=0):
        blocks = []
        lo = min(l for l, r in current_T)
        for b_idx, s in enumerate(starts):
            block_src = current_T[::-1] if (reverse_block_parity and (b_idx % 2 == 1)) else current_T
            shift = s * delta * K - lo
            block = [(l + shift, r + shift) for (l, r) in block_src]
            blocks.append(block)
        return blocks

    # Interleave blocks with parity order control to diversify overlap patterns
    def interleave(blocks, round_idx, do_interleave=True):
        if not do_interleave:
            # sequential (but optionally reverse the block order for added mixing)
            return [iv for blk in blocks for iv in blk]
        maxlen = max(len(b) for b in blocks)
        order = list(range(len(blocks)))
        # Even rounds: forward order; odd rounds: reversed order
        if round_idx % 2 == 1:
            order = order[::-1]
        S = []
        for i in range(maxlen):
            for idx in order:
                blk = blocks[idx]
                if i < len(blk):
                    S.append(blk[i])
        return S

    # Connectors: classic four plus cross3 and cross4 to strengthen long-range coupling.
    def add_connectors(S, starts, delta, K=1):
        s0, s1, s2, s3 = starts
        # scale translations by K to match block placements
        base = delta * K
        connectors = [
            ((s0 - 1) * base, (s1 - 1) * base),   # left cap
            ((s2 + 2) * base, (s3 + 2) * base),   # right cap
            ((s0 + 2) * base, (s2 - 1) * base),   # cross1
            ((s1 + 2) * base, (s3 - 1) * base),   # cross2
            ((s0 + 3) * base, (s3 + 3) * base),   # cross3 (longer)
            ((s0 + 4) * base, (s3 + 4) * base),   # cross4 (new connector)
        ]
        S.extend(connectors)

    # Apply spine rounds
    for ridx in range(depth):
        # pick template
        starts = template_bank[ridx % len(template_bank)] if rotate_starts else template_bank[0]

        # compute current span delta
        lo = min(l for l, r in T)
        hi = max(r for l, r in T)
        delta = hi - lo
        if delta <= 0:
            delta = 1

        # alternate a small spine multiplier K to increase packing density in alternate rounds
        # parity: even rounds get K=2 to create denser translations; odd rounds K=1
        K = 2 if (ridx % 2 == 0) else 1

        # build blocks with optional internal reversal parity
        blocks = build_blocks(T, starts, delta, K=K, reverse_block_parity=reverse_block_parity, round_idx=ridx)

        # decide whether to interleave blocks this round: controlled by interleave_blocks and parity
        do_inter = bool(interleave_blocks and (ridx % 2 == 0))
        S = interleave(blocks, round_idx=ridx, do_interleave=do_inter)

        # add connectors scaled by K
        add_connectors(S, starts, delta, K=K)

        # capacity guard: if expanded S will exceed CAP, break and keep T as-is
        if len(S) > CAP:
            break
        T = S
        if len(T) >= CAP - 16:
            break

    # If already near capacity, return normalized result
    if len(T) >= CAP - 16:
        # normalize to non-negative integers and return
        min_l = min(l for l, _ in T)
        if min_l < 0:
            T = [(l - min_l, r - min_l) for l, r in T]
        return [(int(l), int(r)) for (l, r) in T[:CAP]]

    # Stage 2: micro-phase(s) using delta2 rounds (thin seed sampling)
    def thin_seed(current_T, max_seed):
        if not current_T or max_seed <= 0:
            return []
        n = len(current_T)
        step = max(1, n // max_seed)
        return current_T[::step][:max_seed]

    def build_micro_round(current_T, round_id, budget):
        if not current_T or budget <= 4:
            return []

        glo = min(l for l, r in current_T)
        ghi = max(r for l, r in current_T)
        G = max(1, ghi - glo)

        # delta2 scale: smaller window to avoid blowing global omega
        delta2 = max(1, G // 4)

        # sample thin seed
        per_block_target = max(12, min(64, budget // 6))
        U = thin_seed(current_T, per_block_target)
        if not U:
            return []

        ulo = min(l for l, r in U)

        # four fractional windows inside global span (shifted slightly per round)
        shift_frac = (round_id % 3) * 0.02
        window_fracs = [
            (0.10 + shift_frac, 0.22 + shift_frac),
            (0.33 + shift_frac, 0.45 + shift_frac),
            (0.56 + shift_frac, 0.68 + shift_frac),
            (0.79 + shift_frac, 0.91 + shift_frac),
        ]
        # clamp windows
        window_fracs = [(max(0.05, a), min(0.95, b)) for (a, b) in window_fracs]

        blocks = []
        for idx, (fa, fb) in enumerate(window_fracs):
            win_lo = glo + int(round(fa * G))
            base = win_lo - ulo
            blk = [(l + base, r + base) for (l, r) in U]
            # alternate internal reversal for variety
            if (idx % 2) == (round_id % 2):
                blk = list(reversed(blk))
            blocks.append(blk)

        # interleave micro-blocks; parity alternates order
        block_order = list(range(len(blocks)))
        if round_id % 2 == 1:
            block_order = block_order[::-1]
        micro = []
        maxlen = max(len(b) for b in blocks)
        for i in range(maxlen):
            for idx in block_order:
                blk = blocks[idx]
                if i < len(blk):
                    micro.append(blk[i])

        # add micro-scale connectors including cross4 at delta2 scale
        starts = [2, 6, 10, 14]
        base2 = delta2
        micro_connectors = [
            ((starts[0] - 1) * base2 + glo, (starts[1] - 1) * base2 + glo),
            ((starts[2] + 2) * base2 + glo, (starts[3] + 2) * base2 + glo),
            ((starts[0] + 2) * base2 + glo, (starts[2] - 1) * base2 + glo),
            ((starts[1] + 2) * base2 + glo, (starts[3] - 1) * base2 + glo),
            ((starts[0] + 3) * base2 + glo, (starts[3] + 3) * base2 + glo),  # cross3
            ((starts[0] + 4) * base2 + glo, (starts[3] + 4) * base2 + glo),  # cross4
        ]
        micro.extend(micro_connectors)

        # sprinkle a few long but safe caps (centered) near the micro region to force final FF growth
        mid = glo + G // 2
        cap1 = (mid - max(1, delta2 // 2), mid + max(1, delta2 // 2))
        cap2 = (glo + max(1, delta2 // 3), glo + max(1, int(1.6 * delta2)))
        cap3 = (ghi - max(1, int(1.6 * delta2)), ghi - max(1, delta2 // 3))
        for c in (cap1, cap2, cap3):
            if c[1] > c[0]:
                micro.append(c)

        # enforce budget
        if len(micro) > budget:
            micro = micro[:budget]
        return micro

    steps = min(max(0, int(phase2_iters)), 3)  # allow up to 3 light micro rounds, capacity-guarded
    for mid in range(steps):
        room = CAP - len(T)
        if room <= 8:
            break
        micro = build_micro_round(T, round_id=mid, budget=room)
        if not micro:
            break
        # capacity guard and append
        avail = CAP - len(T)
        if len(micro) > avail:
            micro = micro[:avail]
        # insert micro blocks near tail to interact with late FF colors
        T.extend(micro)

    # Final normalization: shift so minimum left endpoint is >= 0 and convert to integers
    if not T:
        return []
    min_l = min(l for l, r in T)
    if min_l < 0:
        T = [(l - min_l, r - min_l) for l, r in T]

    # Round to integers (safe since constructions use integer math), and trim to CAP
    final = []
    for l, r in T[:CAP]:
        # ensure well-formed integer intervals
        ll = int(round(l))
        rr = int(round(r))
        if rr <= ll:
            rr = ll + 1
        final.append((ll, rr))

    return final

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()