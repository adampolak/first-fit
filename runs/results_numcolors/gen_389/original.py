# EVOLVE-BLOCK-START

def _span_delta(T):
    """Return (lo, hi, delta) for current interval set T."""
    lo = min(l for l, _ in T)
    hi = max(r for _, r in T)
    delta = hi - lo
    if delta <= 0:
        delta = 1
    return lo, hi, delta

def _svc_connectors(starts, delta):
    """Produce four canonical connectors given the start-pattern and block delta."""
    s0, s1, s2, s3 = starts
    return [
        ((s0 - 1) * delta, (s1 - 1) * delta),  # left cap
        ((s2 + 2) * delta, (s3 + 2) * delta),  # right cap
        ((s0 + 2) * delta, (s2 - 1) * delta),  # cross 1
        ((s1 + 2) * delta, (s3 - 1) * delta),  # cross 2
    ]

def _translate_blocks(T, starts, delta, lo):
    """Create four translated blocks from current T and a given start-pattern."""
    blocks = []
    for s in starts:
        base = s * delta - lo
        block = [(l + base, r + base) for (l, r) in T]
        blocks.append(block)
    return blocks

def _assemble_from_blocks(blocks, interleave=False, reverse_order=False, round_seed=0):
    """Flatten blocks into S, optionally interleaving. round_seed adds deterministic rotation."""
    if not interleave:
        if reverse_order:
            blocks = list(reversed(blocks))
        S = []
        for blk in blocks:
            S.extend(blk)
        return S

    # Interleave across blocks with a deterministic per-round rotation
    maxlen = max((len(b) for b in blocks), default=0)
    order = list(range(len(blocks)))
    if reverse_order:
        order = list(reversed(order))
    # rotate order deterministically by round_seed mod len(order)
    if order:
        krot = int(round_seed) % len(order)
        order = order[krot:] + order[:krot]

    S = []
    for i in range(maxlen):
        for idx in order:
            blk = blocks[idx]
            if i < len(blk):
                S.append(blk[i])
    return S

def _apply_round(current_T, starts, do_interleave=False, reverse_order=False, round_seed=0):
    """One KT-round: translate 4 blocks, assemble, then append connectors."""
    lo, hi, delta = _span_delta(current_T)
    # translate 4 blocks
    blocks = _translate_blocks(current_T, starts, delta, lo)
    S = _assemble_from_blocks(blocks, interleave=do_interleave, reverse_order=reverse_order, round_seed=round_seed)
    # append canonical connectors
    S += _svc_connectors(starts, delta)
    return S

def _insert_tail_caps(T, caps, cap_cap):
    """Safely insert tail caps up to a cap_cap budget."""
    room = cap_cap - len(T)
    if room <= 0:
        return T
    for i, iv in enumerate(caps[:room]):
        pos = len(T) - (i * 2 + 1)
        if pos < 0:
            T.append(iv)
        else:
            T.insert(pos, iv)
    return T

def construct_intervals(seed_count=1):
    """
    Deterministic KT-spine with improved deterministic interleaving and a second micro-phase.
    Produces a sequence of intervals for FirstFit with omega kept moderate.
    Returns a list of (l, r) integer pairs with r > l, capped at 9800.
    """
    CAP = 9800
    BASE_SEED = 7  # small deterministic seed used to derive per-round mixing flags

    # Deterministic rotate templates
    template_bank = [
        (2, 6, 10, 14),  # classic KT
        (1, 5, 9, 13),   # left-shifted
        (3, 7, 11, 15),  # right-shifted
        (4, 8, 12, 16),  # stretched-right
    ]

    # Seed sequence
    if seed_count <= 1:
        T = [(0, 1)]
    else:
        step = 3
        seeds = min(4, max(1, int(seed_count)))
        T = [(i * step, i * step + 1) for i in range(seeds)]

    # Stage 1: six KT-style rounds with deterministic per-round interleaving derived from BASE_SEED
    for ridx in range(6):
        starts = template_bank[ridx % len(template_bank)]
        nxt_size = 4 * len(T) + 4
        if nxt_size > CAP:
            break
        # derive deterministic per-round flags from BASE_SEED and round index
        round_val = (BASE_SEED * 31 + ridx * 17) & 0xFFFF
        # enforce interleaving every round to maximize cross-block mixing; keep deterministic reversal
        do_interleave = True
        reverse_order = ((round_val // 2) % 2 == 1)
        # supply round_seed to _apply_round to rotate interleaving order deterministically
        T = _apply_round(T, starts, do_interleave=do_interleave, reverse_order=reverse_order, round_seed=round_val)
        if len(T) >= CAP:
            return T[:CAP]

    # Stage 1b (optional): early return guard
    if len(T) >= CAP - 16:
        lo, hi, delta = _span_delta(T)
        long_connectors = [
            (lo + max(1, int(round(0.05 * delta))), lo + max(1, int(round(0.95 * delta)))),
            (lo + max(1, int(round(0.15 * delta))), lo + max(1, int(round(0.85 * delta)))),
            (lo + max(1, int(round(0.30 * delta))), lo + max(1, int(round(0.70 * delta)))),
        ]
        for c in long_connectors:
            if len(T) < CAP:
                T.append(c)
        return T[:CAP]

    # Micro-phase A: tail caps to mix colors without exploding omega
    lo, hi, delta = _span_delta(T)
    tail_caps = [
        (lo + max(1, int(round(delta * 0.08))), lo + max(1, int(round(delta * 0.60)))),
        (lo + max(1, int(round(delta * 0.25))), lo + max(1, int(round(delta * 0.75)))),
        (lo + max(1, int(round(delta * 0.75))), lo + max(1, int(round(delta * 0.92)))),
    ]
    T = _insert_tail_caps(T, tail_caps, CAP)

    if len(T) >= CAP - 16:
        return T[:CAP]

    # Micro-phase B: delta2 micro-rounds with thin seeds
    def _thin_seed(current_T, max_seed):
        n = len(current_T)
        if n == 0 or max_seed <= 0:
            return []
        step = max(1, n // max_seed)
        return current_T[::step][:max_seed]

    def _build_micro_round(current_T, budget, iter_id=0, alt=False, round_seed=0):
        if not current_T or budget <= 8:
            return []

        glo = min(l for l, _ in current_T)
        ghi = max(r for _, r in current_T)
        G = max(1, ghi - glo)

        # seed
        seed_sz = max(8, min(40, len(current_T) // 250))
        if alt:
            seed_sz = max(8, min(64, len(current_T) // 200))
        U = _thin_seed(current_T, seed_sz)
        if not U:
            return []

        ulo = min(l for l, _ in U)

        if not alt:
            shift = ((round_seed + iter_id) % 3) * 0.02
            window_fracs = [
                (0.12 + shift, 0.22 + shift),
                (0.35 + shift, 0.45 + shift),
                (0.58 + shift, 0.68 + shift),
                (0.80 + shift, 0.90 + shift),
            ]
            window_fracs = [
                (max(0.05, min(0.90, a)),
                 max(0.10, min(0.95, b)))
                for (a, b) in window_fracs
            ]
        else:
            window_fracs = [(0.05, 0.15), (0.28, 0.38), (0.60, 0.70), (0.82, 0.92)]

        blocks = []
        for wi, (fa, fb) in enumerate(window_fracs):
            win_lo = glo + int(round(fa * G))
            base = win_lo - ulo
            block = [(l + base, r + base) for (l, r) in U]
            tag = iter_id + (round_seed % 5)
            if ((int(round(fa * 100)) // 5) + tag + wi) % 2 == 1:
                block = list(reversed(block))
            blocks.append(block)

        micro = _assemble_from_blocks(blocks, interleave=True, reverse_order=(iter_id % 2 == 1), round_seed=round_seed)

        micro_connectors = [
            (glo + int(round(0.08 * G)), glo + int(round(0.30 * G))),
            (glo + int(round(0.60 * G)), glo + int(round(0.92 * G))),
            (glo + int(round(0.26 * G)), glo + int(round(0.56 * G))),
            (glo + int(round(0.44 * G)), glo + int(round(0.78 * G))),
        ]
        if alt:
            micro_connectors.append((glo + int(round(0.18 * G)), glo + int(round(0.84 * G))))
            micro_connectors.append((glo + int(round(0.10 * G)), glo + int(round(0.88 * G))))
        for a, b in micro_connectors:
            if b > a:
                micro.append((a, b))

        if len(micro) > budget:
            micro = micro[:budget]
        return micro

    # Run two main micro rounds (safe cap), with deterministic seeds
    steps = 2
    for it in range(steps):
        room = CAP - len(T)
        if room <= 8:
            break
        round_val = (BASE_SEED * 97 + it * 43) & 0xFFFF
        micro = _build_micro_round(T, room, iter_id=it, alt=False, round_seed=round_val)
        if not micro:
            break
        if len(micro) > room:
            micro = micro[:room]
        T.extend(micro)

    # Alternate micro family
    room = CAP - len(T)
    if room > 8:
        round_val = (BASE_SEED * 113 + 2 * 59) & 0xFFFF
        micro_alt = _build_micro_round(T, room, iter_id=2, alt=True, round_seed=round_val)
        if micro_alt:
            if len(micro_alt) > room:
                micro_alt = micro_alt[:room]
            T.extend(micro_alt)

    # Stage 3: small pinning micro-phase (short intervals to cap omega)
    if len(T) < CAP - 8:
        glo2 = min(l for l, _ in T)
        ghi2 = max(r for _, r in T)
        G2 = max(1, ghi2 - glo2)
        U2 = _thin_seed(T, max(8, min(24, len(T) // 400)))
        if U2:
            ulo2 = min(l for l, r in U2)
            eps = max(1, G2 // 512)
            windows2 = [(0.05, 0.15), (0.28, 0.38), (0.60, 0.70), (0.82, 0.92)]
            blocks2 = []
            for (fa, fb) in windows2:
                win_lo2 = glo2 + int(round(fa * G2))
                base2 = win_lo2 - ulo2
                block = []
                for idx, (l, r) in enumerate(U2):
                    mid = (l + r) // 2
                    L = mid + base2 - (eps // 2) + (idx % 3)
                    R = L + eps
                    if R > L:
                        block.append((L, R))
                blocks2.append(block)

            micro2 = []
            if blocks2:
                maxlen2 = max(len(b) for b in blocks2)
                order2 = [3, 1, 0, 2]
                for i in range(maxlen2):
                    for idx in order2:
                        blk = blocks2[idx]
                        if i < len(blk):
                            micro2.append(blk[i])

            caps2 = [
                (glo2 + int(round(0.09 * G2)), glo2 + int(round(0.09 * G2)) + eps),
                (glo2 + int(round(0.64 * G2)), glo2 + int(round(0.64 * G2)) + eps),
            ]
            for a, b in caps2:
                if b > a:
                    micro2.append((a, b))

            room2 = CAP - len(T)
            if room2 > 0:
                if len(micro2) > room2:
                    micro2 = micro2[:room2]
                T.extend(micro2)

    # Stage 4: tower-and-cap micro-phase (five towers with caps)
    if len(T) < CAP - 12:
        glo3 = min(l for l, r in T)
        ghi3 = max(r for l, r in T)
        G3 = max(1, ghi3 - glo3)

        U3 = _thin_seed(T, max(12, min(36, len(T) // 300)))
        if U3:
            ulo3 = min(l for l, r in U3)
            eps3 = max(1, G3 // 768)

            towers = [(0.10, 0.18), (0.28, 0.36), (0.46, 0.54), (0.64, 0.72), (0.82, 0.90)]
            levels = 3
            tower_blocks = []
            for t_idx, (fa, fb) in enumerate(towers):
                win_lo3 = glo3 + int(round(fa * G3))
                base3 = win_lo3 - ulo3
                block = []
                for idx, (l, r) in enumerate(U3):
                    mid = (l + r) // 2
                    for lv in range(levels):
                        jitter = (idx + lv + t_idx) % 4
                        L = mid + base3 - (eps3 // 2) + jitter + lv
                        R = L + eps3
                        if R > L:
                            block.append((L, R))
                tower_blocks.append(block)

            tower_phase = []
            if tower_blocks:
                maxlen3 = max(len(b) for b in tower_blocks)
                order3 = [0, 2, 4, 1, 3]
                for i in range(maxlen3):
                    for idx in order3:
                        blk = tower_blocks[idx]
                        if i < len(blk):
                            tower_phase.append(blk[i])

            caps3 = [
                (glo3 + int(round(0.15 * G3)), glo3 + int(round(0.45 * G3))),
                (glo3 + int(round(0.38 * G3)), glo3 + int(round(0.70 * G3))),
                (glo3 + int(round(0.60 * G3)), glo3 + int(round(0.88 * G3))),
            ]
            for a, b in caps3:
                if b > a:
                    tower_phase.append((a, b))

            room3 = CAP - len(T)
            if room3 > 0 and tower_phase:
                if len(tower_phase) > room3:
                    tower_phase = tower_phase[:room3]
                T.extend(tower_phase)

    # SECOND MICRO-PHASE (compact, deterministic, cross-scale connectors)
    room = CAP - len(T)
    if room > 8:
        glo4 = min(l for l, _ in T)
        ghi4 = max(r for _, r in T)
        G4 = max(1, ghi4 - glo4)
        window_fracs_b = [(0.06, 0.16), (0.30, 0.40), (0.54, 0.64), (0.78, 0.88)]
        U4 = _thin_seed(T, max(8, min(28, len(T) // 260)))
        if U4:
            ulo4 = min(l for l, _ in U4)
            blocks_b = []
            for wi, (fa, fb) in enumerate(window_fracs_b):
                win_lo = glo4 + int(round(fa * G4))
                base = win_lo - ulo4
                block = [(l + base, r + base) for (l, r) in U4]
                if (BASE_SEED + wi) % 2 == 0:
                    block = list(reversed(block))
                blocks_b.append(block)
            micro_b = _assemble_from_blocks(blocks_b, interleave=True, reverse_order=False, round_seed=(BASE_SEED * 3))
            # add a few conservative long connectors tuned to full span to couple colors across the backbone
            cross_connectors = [
                (glo4 + int(round(0.12 * G4)), glo4 + int(round(0.66 * G4))),
                (glo4 + int(round(0.20 * G4)), glo4 + int(round(0.82 * G4))),
                (glo4 + int(round(0.36 * G4)), glo4 + int(round(0.90 * G4))),
            ]
            for a, b in cross_connectors:
                if b > a:
                    micro_b.append((a, b))
            if len(micro_b) > room:
                micro_b = micro_b[:room]
            T.extend(micro_b)

    # Final capacity guard
    if len(T) > CAP:
        T = T[:CAP]

    # Normalize endpoints to integers with r > l
    if not T:
        return []
    min_l = min(l for l, _ in T)
    if min_l < 0:
        T = [(l - min_l, r - min_l) for (l, r) in T]

    intervals = []
    for (l, r) in T:
        li = int(l)
        ri = int(r)
        if ri <= li:
            ri = li + 1
        intervals.append((li, ri))

    if len(intervals) > CAP:
        intervals = intervals[:CAP]
    return intervals

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()