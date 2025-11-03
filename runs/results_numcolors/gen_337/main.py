# EVOLVE-BLOCK-START

def construct_intervals(seed_count=1):
    """
    Build an interval sequence to push FirstFit colors while keeping omega near 10.
    Inputs/outputs preserved relative to prior versions.

    Returns:
      list of (l, r) integer pairs with r > l (open intervals), length <= 9800.
    """
    # -------------------
    # Configuration knobs
    # -------------------
    CAP = 9800
    ROUNDS = 6
    BASE_SEED = 0x5a77b  # fixed deterministic seed

    # Six-template spine bank (expanded over earlier 4 pattern variants)
    TEMPLATE_BANK = [
        (2, 6, 10, 14),  # classic KT
        (1, 5, 9, 13),   # left-shifted
        (3, 7, 11, 15),  # right-shifted
        (4, 8, 12, 16),  # stretched-right
        (2, 4, 8, 12),   # compressed-left
        (3, 5, 9, 13),   # gentle-left
    ]

    # Micro window families (two-phase)
    WINDOWS_A = [(0.12, 0.22), (0.35, 0.45), (0.58, 0.68), (0.80, 0.90)]
    WINDOWS_B = [(0.06, 0.16), (0.30, 0.40), (0.54, 0.64), (0.78, 0.88)]

    # -----------------
    # Utility functions
    # -----------------
    def _span_delta(T):
        lo = min(l for l, _ in T)
        hi = max(r for _, r in T)
        delta = hi - lo
        if delta <= 0:
            delta = 1
        return lo, hi, delta

    def _det_hash(x):
        """Deterministic integer hash."""
        x ^= (x >> 17)
        x *= 0xed5ad4bb
        x &= 0xffffffff
        x ^= (x >> 11)
        x *= 0xac4c1b51
        x &= 0xffffffff
        x ^= (x >> 15)
        x *= 0x31848bab
        x &= 0xffffffff
        x ^= (x >> 14)
        return x & 0xffffffff

    def _choose_template(templates, base_seed, round_idx):
        idx = _det_hash(base_seed ^ (round_idx * 0x9e3779b1)) % len(templates)
        return templates[idx]

    def _translate_block(T, shift):
        return [(l + shift, r + shift) for (l, r) in T]

    def _assemble_blocks(blocks, interleave=False, reverse=False, rotation=0):
        if not blocks:
            return []
        if not interleave:
            use = list(reversed(blocks)) if reverse else blocks
            out = []
            for b in use:
                out.extend(b)
            return out
        # interleave with deterministic block order rotation
        order = list(range(len(blocks)))
        if reverse:
            order.reverse()
        if order:
            k = rotation % len(order)
            order = order[k:] + order[:k]
        out = []
        maxlen = max(len(b) for b in blocks)
        for i in range(maxlen):
            for idx in order:
                blk = blocks[idx]
                if i < len(blk):
                    out.append(blk[i])
        return out

    def _classic_connectors(starts, delta):
        s0, s1, s2, s3 = starts
        return [
            ((s0 - 1) * delta, (s1 - 1) * delta),  # left cap
            ((s2 + 2) * delta, (s3 + 2) * delta),  # right cap
            ((s0 + 2) * delta, (s2 - 1) * delta),  # cross 1
            ((s1 + 2) * delta, (s3 - 1) * delta),  # cross 2
        ]

    def _long_range_connectors(lo, hi, base_seed, count=5):
        """Deterministic long-range connectors anchored by fractional span."""
        span = max(1, hi - lo)
        # Fractions engineered to cover left/mid/right with overlaps but not all at once
        base_fracs = [
            (0.10, 0.66),
            (0.18, 0.84),
            (0.32, 0.90),
            (0.44, 0.78),
            (0.22, 0.72),
            (0.12, 0.56),
        ]
        out = []
        for i, (a, b) in enumerate(base_fracs[:max(0, count)]):
            # Small deterministic per-connector jitter to diversify spans
            jitter = ((_det_hash(base_seed + i) % 5) - 2) * 0.005  # in [-0.01, 0.01]
            A = max(0.02, min(0.92, a + jitter))
            B = max(0.08, min(0.96, b + jitter))
            L = lo + max(1, int(round(A * span)))
            R = lo + max(1, int(round(B * span)))
            if R <= L:
                R = L + 1
            out.append((L, R))
        return out

    def _thin_even_seed(T, max_seed, lb=8, ub=64):
        target = max(lb, min(ub, max_seed))
        n = len(T)
        if n == 0:
            return []
        step = max(1, n // target)
        seed = T[::step][:target]
        if not seed:
            seed = list(T)
        return seed

    def _micro_round(current_T, window_fracs, base_seed, iter_tag, budget, add_bridge=False):
        """Fractional-window micro round with deterministic interleaving and connectors."""
        if not current_T or budget <= 8:
            return []
        glo = min(l for l, _ in current_T)
        ghi = max(r for _, r in current_T)
        G = max(1, ghi - glo)

        seed_sz = max(8, min(48, len(current_T) // 240))
        U = _thin_even_seed(current_T, seed_sz, lb=8, ub=48)
        if not U:
            return []
        ulo = min(l for l, _ in U)

        # Deterministic window shift for diversity
        shift_unit = (_det_hash(base_seed ^ (iter_tag * 0x632b)) % 7) * 0.005  # 0, .005, .010, ..., .030
        blocks = []
        for widx, (fa, fb) in enumerate(window_fracs):
            a = max(0.05, min(0.90, fa + shift_unit))
            win_lo = glo + int(round(a * G))
            base = win_lo - ulo
            blk = _translate_block(U, base)
            # Deterministic internal reversal to break symmetry
            if ((_det_hash(base_seed + widx + iter_tag) >> 3) & 1) == 1:
                blk = list(reversed(blk))
            blocks.append(blk)

        # Interleave with deterministic order (even iter -> forward; odd -> reversed)
        rotation = (_det_hash(base_seed + iter_tag) % 4)
        micro = _assemble_blocks(
            blocks,
            interleave=True,
            reverse=bool(iter_tag % 2),
            rotation=rotation
        )

        # Fractional connectors at micro scale
        micro_connectors = [
            (glo + int(round(0.08 * G)), glo + int(round(0.30 * G))),
            (glo + int(round(0.60 * G)), glo + int(round(0.92 * G))),
            (glo + int(round(0.26 * G)), glo + int(round(0.56 * G))),
            (glo + int(round(0.44 * G)), glo + int(round(0.78 * G))),
        ]
        if add_bridge:
            micro_connectors.append((glo + int(round(0.18 * G)), glo + int(round(0.84 * G))))
        for a, b in micro_connectors:
            if b > a:
                micro.append((a, b))

        if len(micro) > budget:
            micro = micro[:budget]
        return micro

    def _spine_density_boost(T, base_seed, cap_left):
        """
        Small CAP-aware density booster that injects short pins in safe windows near edges.
        Keeps omega modest by using very short intervals and non-overlapping windows.
        """
        if cap_left <= 16 or not T:
            return []
        glo, ghi, G = _span_delta(T)
        eps = max(1, G // 1024)
        # Windows near edges (sparse), chosen to avoid central high-density regions
        edge_windows = [(0.06, 0.10), (0.90, 0.94)]
        U = _thin_even_seed(T, max_seed=max(12, min(40, len(T) // 300)), lb=8, ub=40)
        pins = []
        for widx, (fa, fb) in enumerate(edge_windows):
            anchor = glo + int(round(((fa + fb) / 2.0) * G))
            jitter_base = (_det_hash(base_seed + widx) % 5) - 2  # -2..+2
            for sidx, (l, r) in enumerate(U):
                mid = (l + r) // 2
                shift = (sidx % 7) + jitter_base
                L = anchor + (mid % 5) + shift
                R = L + eps
                if R > L:
                    pins.append((L, R))
                if len(pins) >= cap_left // 8:  # keep this boost tiny per round
                    break
            if len(pins) >= cap_left // 8:
                break
        if len(pins) > cap_left:
            pins = pins[:cap_left]
        return pins

    # -------------------
    # Seed initialization
    # -------------------
    if seed_count <= 1:
        T = [(0, 1)]
    else:
        # Conservative multi-seed: staggered to avoid early omega spikes
        seeds = min(4, max(1, int(seed_count)))
        step = 3
        T = [(i * step, i * step + 1) for i in range(seeds)]

    # -------------------------------
    # Stage 1: Hex-template KT spine
    # -------------------------------
    for ridx in range(ROUNDS):
        # Capacity pre-check for KT growth: size -> 4*size + 4
        if 4 * len(T) + 4 > CAP:
            break

        starts = _choose_template(TEMPLATE_BANK, BASE_SEED, ridx)
        lo, hi, delta = _span_delta(T)

        # Build 4 blocks via translation
        blocks = []
        for s in starts:
            shift = s * delta - lo
            blocks.append(_translate_block(T, shift))

        # Controlled interleaving: even -> interleave with rotation; odd -> reversed sequential
        interleave = (ridx % 2 == 0)
        reverse = (ridx % 2 == 1)
        rotation = _det_hash(BASE_SEED ^ (ridx * 0x5851f42d)) % 4
        S = _assemble_blocks(blocks, interleave=interleave, reverse=reverse, rotation=rotation)

        # Append classic SVC connectors
        S.extend(_classic_connectors(starts, delta))

        # CAP-aware density boost per round (tiny)
        cap_left = CAP - len(S)
        if cap_left > 16:
            pins = _spine_density_boost(S, BASE_SEED + ridx * 17, cap_left)
            if pins:
                S.extend(pins)

        T = S
        if len(T) >= CAP:
            T = T[:CAP]
            return _normalize_intervals(T, CAP)

    # -------------------------------
    # Post-spine long-range connectors
    # -------------------------------
    if len(T) < CAP - 8:
        lo, hi, _ = _span_delta(T)
        # 4–6 deterministic long connectors (gated by CAP)
        long_cons = _long_range_connectors(lo, hi, BASE_SEED, count=6)
        room = CAP - len(T)
        if long_cons:
            T.extend(long_cons[:room])

    if len(T) >= CAP - 8:
        return _normalize_intervals(T[:CAP], CAP)

    # -----------------------------------------------
    # CAP-aware two-phase micro-phase budget planning
    # -----------------------------------------------
    remaining = CAP - len(T)
    if remaining <= 8:
        return _normalize_intervals(T, CAP)

    # Reserve two micro passes (A then B); split deterministically 60/40
    budget_A = int(remaining * 0.60)
    budget_B = remaining - budget_A
    budget_A = max(0, budget_A)
    budget_B = max(0, budget_B)

    # ------------------------
    # Micro-phase A (WINDOWS_A)
    # ------------------------
    if budget_A > 8:
        microA = _micro_round(T, WINDOWS_A, BASE_SEED ^ 0x1111, iter_tag=0, budget=budget_A, add_bridge=False)
        if microA:
            T.extend(microA[:budget_A])

    # Capacity refresh
    remaining = CAP - len(T)
    if remaining <= 8:
        return _normalize_intervals(T, CAP)

    # ------------------------
    # Micro-phase B (WINDOWS_B)
    # ------------------------
    if budget_B > 8:
        microB = _micro_round(T, WINDOWS_B, BASE_SEED ^ 0x2222, iter_tag=1, budget=budget_B, add_bridge=True)
        if microB:
            T.extend(microB[:budget_B])

    # Final clamp and normalize
    if len(T) > CAP:
        T = T[:CAP]
    return _normalize_intervals(T, CAP)


# ---- Normalization helper kept outside construct_intervals for clarity ----
def _normalize_intervals(T, CAP):
    """Normalize to integer endpoints with r > l and clamp length to CAP."""
    if not T:
        return []
    # Translate if needed to keep non-negative
    min_l = min(l for l, _ in T)
    if min_l < 0:
        T = [(l - min_l, r - min_l) for (l, r) in T]
    out = []
    for (l, r) in T:
        li = int(l)
        ri = int(r)
        if ri <= li:
            ri = li + 1
        out.append((li, ri))
    if len(out) > CAP:
        out = out[:CAP]
    return out

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()