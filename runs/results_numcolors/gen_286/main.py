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

def _assemble_from_blocks(blocks, interleave=False, reverse_order=False):
    """Flatten blocks into S, optionally interleaving."""
    if not interleave:
        if reverse_order:
            blocks = list(reversed(blocks))
        S = []
        for blk in blocks:
            S.extend(blk)
        return S

    # Interleave across blocks
    maxlen = max((len(b) for b in blocks), default=0)
    order = list(range(len(blocks)))
    if reverse_order:
        order = list(reversed(order))
    S = []
    for i in range(maxlen):
        for idx in order:
            blk = blocks[idx]
            if i < len(blk):
                S.append(blk[i])
    return S

def _apply_round(current_T, starts, do_interleave=False, reverse_order=False, pre_caps=False, zigzag_inside=False):
    """One KT-round: translate 4 blocks, optionally pre-cap, assemble, then append remaining connectors."""
    lo, hi, delta = _span_delta(current_T)
    # translate 4 blocks
    blocks = _translate_blocks(current_T, starts, delta, lo)

    # optional zigzag inside-block ordering to disturb FF reuse
    if zigzag_inside and blocks:
      blocks = [blk if (i % 2 == 0) else list(reversed(blk)) for i, blk in enumerate(blocks)]

    # assemble blocks
    S = _assemble_from_blocks(blocks, interleave=do_interleave, reverse_order=reverse_order)

    # connectors with optional CAP-first gating: place caps before blocks, crosses after
    conns = _svc_connectors(starts, delta)
    if pre_caps:
        # left and right caps first; cross connectors after blocks
        left_cap, right_cap, cross1, cross2 = conns
        S = [left_cap, right_cap] + S + [cross1, cross2]
    else:
        # classic post-append of all four connectors
        S += conns
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

def _add_belt_connectors(T, cap_cap):
    """
    Insert a small belt of long-range cross-scale connectors near the tail.
    Fractions chosen to avoid a dense single-point core while coupling distant regions.
    """
    if not T:
        return T
    room = cap_cap - len(T)
    if room <= 0:
        return T
    glo, ghi, G = _span_delta(T)
    # Deterministic fractional-span connectors (belt)
    fracs = [(0.06, 0.34), (0.28, 0.64), (0.14, 0.54), (0.54, 0.88)]
    belt = []
    for a, b in fracs:
        L = glo + int(round(max(0.0, min(0.95, a)) * G))
        R = glo + int(round(max(0.05, min(0.98, b)) * G))
        if R > L:
            belt.append((L, R))
    # Insert near the tail to intersect many active colors
    return _insert_tail_caps(T, belt, cap_cap)

def construct_intervals(seed_count=1):
    """
    Deterministic KT-spine with modular, readable structure.
    Produces a sequence of intervals for FirstFit with omega kept moderate.
    Returns a list of (l, r) integer pairs with r > l, capped at 9800.
    """
    CAP = 9800

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

    # Stage 1: six KT-style rounds with alternating interleave/reverse order
    for ridx in range(6):
        starts = template_bank[ridx % len(template_bank)]
        nxt_size = 4 * len(T) + 4
        if nxt_size > CAP:
            break
        do_interleave = (ridx % 2 == 0)      # even rounds interleave
        reverse_order = (ridx % 2 == 1)      # odd rounds reverse order
        pre_caps = (ridx % 2 == 0)           # CAP-first gating on even rounds
        zigzag_inside = (ridx % 2 == 1)      # zigzag inside blocks on odd rounds
        T = _apply_round(T, starts, do_interleave=do_interleave, reverse_order=reverse_order,
                         pre_caps=pre_caps, zigzag_inside=zigzag_inside)
        if len(T) >= CAP:
            return T[:CAP]

    # Insert a small belt of cross-scale connectors before micro-phases
    T = _add_belt_connectors(T, CAP)

    # Stage 1b (optional): early return guard
    if len(T) >= CAP - 16:
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

    # Stage 2: delta2 micro-rounds with thin seeds
    def _thin_seed(current_T, max_seed):
        n = len(current_T)
        if n == 0 or max_seed <= 0:
            return []
        step = max(1, n // max_seed)
        return current_T[::step][:max_seed]

    def _build_micro_round(current_T, budget, iter_id=0, alt=False):
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
            shift = (iter_id % 3) * 0.02
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
        for (fa, fb) in window_fracs:
            win_lo = glo + int(round(fa * G))
            base = win_lo - ulo
            block = [(l + base, r + base) for (l, r) in U]
            # small parity perturbation
            tag = iter_id if not alt else (iter_id + 1)
            if ((int(round(fa * 100)) // 5) + tag) % 2 == 1:
                block = list(reversed(block))
            blocks.append(block)

        micro = _assemble_from_blocks(blocks, interleave=True, reverse_order=(iter_id % 2 == 1))

        # connectors
        micro_connectors = [
            (glo + int(round(0.08 * G)), glo + int(round(0.30 * G))),
            (glo + int(round(0.60 * G)), glo + int(round(0.92 * G))),
            (glo + int(round(0.26 * G)), glo + int(round(0.56 * G))),
            (glo + int(round(0.44 * G)), glo + int(round(0.78 * G))),
        ]
        for a, b in micro_connectors:
            if b > a:
                micro.append((a, b))

        if len(micro) > budget:
            micro = micro[:budget]
        return micro

    # Run two main micro rounds (safe cap)
    steps = 2
    for it in range(steps):
        room = CAP - len(T)
        if room <= 8:
            break
        micro = _build_micro_round(T, room, iter_id=it, alt=False)
        if not micro:
            break
        if len(micro) > room:
            micro = micro[:room]
        T.extend(micro)

    # Alternate micro family
    room = CAP - len(T)
    if room > 8:
        micro_alt = _build_micro_round(T, room, iter_id=2, alt=True)
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