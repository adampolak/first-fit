# EVOLVE-BLOCK-START

from typing import List, Tuple

def _span_delta(T: List[Tuple[int, int]]) -> Tuple[int, int, int]:
    """Return (lo, hi, delta) for current integer interval set T."""
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

def _assemble_braided(blocks, round_seed: int) -> List[Tuple[int, int]]:
    """
    Braided assembly: interleave using stripes of length k that rotates per round.
    This differs from simple interleave by using varying stripe widths and rotated order.
    """
    if not blocks:
        return []
    m = max(len(b) for b in blocks)
    order = list(range(len(blocks)))
    # Deterministic rotation of block order by round seed
    krot = (round_seed // 7) % len(order)
    order = order[krot:] + order[:krot]

    # Stripe size rotates in {1,2,3,4} deterministically per round
    stripe = 1 + (round_seed % 4)

    S = []
    # Interleave with stripes: take 'stripe' intervals from each block before moving to next
    idxs = [0, 0, 0, 0]
    remaining = True
    while remaining:
        remaining = False
        for bi in order:
            blk = blocks[bi]
            i = idxs[bi]
            if i < len(blk):
                remaining = True
                j = min(i + stripe, len(blk))
                S.extend(blk[i:j])
                idxs[bi] = j
    return S

def _apply_braided_round(current_T, starts, round_seed: int):
    """One braided round: translate 4 blocks, assemble with braided mixer, then append connectors."""
    lo, hi, delta = _span_delta(current_T)
    blocks = _translate_blocks(current_T, starts, delta, lo)
    S = _assemble_braided(blocks, round_seed=round_seed)
    S += _svc_connectors(starts, delta)
    return S

def _thin_seed(current_T, max_seed):
    """Evenly spaced seed sample from current_T."""
    n = len(current_T)
    if n == 0 or max_seed <= 0:
        return []
    step = max(1, n // max_seed)
    return current_T[::step][:max_seed]

def _densify_after_round(T, cap, base_seed: int, round_idx: int) -> List[Tuple[int, int]]:
    """
    CAP-bounded density booster after a spine round.
    Inserts small translated seeds into 2 windows per round, deterministically chosen.
    """
    room = cap - len(T)
    if room <= 12:
        return T

    lo, hi, delta = _span_delta(T)
    G = max(1, hi - lo)

    # Choose small seed; keep thin to avoid blowing up omega
    seed_sz = max(8, min(28, len(T) // 350))
    U = _thin_seed(T, seed_sz)
    if not U:
        return T
    ulo = min(l for l, _ in U)

    # Deterministic window selection per round
    # Two windows per round, rotated across a bank of 6
    bank = [
        (0.14, 0.22), (0.34, 0.42),
        (0.56, 0.64), (0.76, 0.84),
        (0.20, 0.28), (0.68, 0.76)
    ]
    w1 = bank[(base_seed + 3 * round_idx) % len(bank)]
    w2 = bank[(base_seed + 5 * round_idx + 2) % len(bank)]
    windows = [w1, w2]

    # Budget per densifier pass: small and CAP-bounded
    budget = min(room, seed_sz * len(windows))
    micro = []
    for wi, (fa, fb) in enumerate(windows):
        win_lo = lo + int(round(fa * G))
        base = win_lo - ulo
        block = [(l + base, r + base) for (l, r) in U]
        # reverse every other window for symmetry break
        if ((base_seed + round_idx + wi) % 2) == 1:
            block = list(reversed(block))
        micro.extend(block)

    if len(micro) > budget:
        micro = micro[:budget]
    return T + micro

def _long_connectors(T, cap, count=5):
    """Add a small deterministic set of long-range connectors to couple colors across the backbone."""
    room = cap - len(T)
    if room <= 0:
        return T
    lo, hi, _ = _span_delta(T)
    G = max(1, hi - lo)
    # A short bank of cross-intervals spanning long ranges, positioned safely inside the span
    fracs = [
        (0.10, 0.62),
        (0.18, 0.78),
        (0.32, 0.86),
        (0.44, 0.90),
        (0.22, 0.70),
        (0.58, 0.92)
    ][:count]
    connectors = []
    for a, b in fracs:
        L = lo + int(round(a * G))
        R = lo + int(round(b * G))
        if R <= L:
            R = L + 1
        connectors.append((L, R))
    add = connectors[:room]
    return T + add

def _micro_phase(current_T, cap, base_seed: int, windows, seed_factor=300, max_seed=40):
    """Generic micro-phase builder with given window fractions."""
    room = cap - len(current_T)
    if room <= 8:
        return current_T

    glo, ghi, _ = _span_delta(current_T)
    G = max(1, ghi - glo)

    # Thin, evenly-spaced seed; factor guards micro density
    seed_sz = max(8, min(max_seed, len(current_T) // max(1, seed_factor)))
    U = _thin_seed(current_T, seed_sz)
    if not U:
        return current_T
    ulo = min(l for l, _ in U)

    blocks = []
    for wi, (fa, fb) in enumerate(windows):
        win_lo = glo + int(round(fa * G))
        base = win_lo - ulo
        block = [(l + base, r + base) for (l, r) in U]
        # symmetric breaking via deterministic reversal
        if ((base_seed + wi) % 2) == 0:
            block = list(reversed(block))
        blocks.append(block)

    # Interleave the micro blocks with a fixed pattern that differs from the spine braiding
    # Use order [0,2,1,3] for mixing
    order = [0, 2, 1, 3][:len(blocks)]
    micro = []
    maxlen = max(len(b) for b in blocks)
    for i in range(maxlen):
        for idx in order:
            blk = blocks[idx]
            if i < len(blk):
                micro.append(blk[i])

    # Add micro-scale connectors (fractional-span analog of KT caps), conservative to guard omega
    micro_connectors = [
        (glo + int(round(0.08 * G)), glo + int(round(0.30 * G))),
        (glo + int(round(0.64 * G)), glo + int(round(0.92 * G))),
        (glo + int(round(0.26 * G)), glo + int(round(0.56 * G))),
        (glo + int(round(0.44 * G)), glo + int(round(0.78 * G))),
    ]
    for a, b in micro_connectors:
        if b > a:
            micro.append((a, b))

    if len(micro) > room:
        micro = micro[:room]
    return current_T + micro

def construct_intervals(seed_count=1):
    """
    Braided six-template backbone with mandatory interleaving, density boosters,
    long connectors, and two micro-phase passes. Returns a list of (l, r) integer pairs.
    """
    CAP = 9800
    BASE_SEED = 173  # seed for deterministic rotations

    # Six strong start-pattern templates, rotated deterministically
    template_bank = [
        (2, 6, 10, 14),  # classic KT
        (1, 5, 9, 13),   # left-shifted
        (3, 7, 11, 15),  # right-shifted
        (4, 8, 12, 16),  # stretched-right
        (2, 4, 8, 12),   # compressed-left
        (3, 5, 9, 13),   # gentle-left
    ]

    # Seed: a single unit interval keeps omega minimal in early growth
    T = [(0, 1)] if seed_count <= 1 else [(i * 3, i * 3 + 1) for i in range(min(4, max(1, int(seed_count))))]

    # Stage 1: six braided rounds with deterministic rotation and mandatory interleaving
    for ridx in range(6):
        # Pre-check capacity: allow the core round to fit
        nxt_size = 4 * len(T) + 4
        if nxt_size > CAP:
            break
        starts = template_bank[(BASE_SEED + ridx) % len(template_bank)]
        round_seed = (BASE_SEED * 97 + ridx * 131) & 0x7fffffff
        T = _apply_braided_round(T, starts, round_seed=round_seed)
        if len(T) >= CAP:
            T = T[:CAP]
            return T
        # CAP-bounded densification after each round (small boost)
        T = _densify_after_round(T, CAP, base_seed=BASE_SEED, round_idx=ridx)
        if len(T) >= CAP:
            T = T[:CAP]
            return T

    # Deterministic long-range connectors to couple distant colors
    T = _long_connectors(T, CAP, count=6)
    if len(T) >= CAP - 8:
        return T[:CAP]

    # Stage 2A: micro-phase with window family A
    windows_A = [(0.12, 0.22), (0.35, 0.45), (0.58, 0.68), (0.80, 0.90)]
    T = _micro_phase(T, CAP, base_seed=BASE_SEED + 11, windows=windows_A, seed_factor=280, max_seed=44)
    if len(T) >= CAP - 8:
        return T[:CAP]

    # Stage 2B: second independent micro-phase with distinct windows B
    windows_B = [(0.06, 0.16), (0.30, 0.40), (0.54, 0.64), (0.78, 0.88)]
    T = _micro_phase(T, CAP, base_seed=BASE_SEED + 37, windows=windows_B, seed_factor=260, max_seed=48)
    if len(T) > CAP:
        T = T[:CAP]

    # Normalize endpoints to integers with r > l and non-negative coordinates
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
        if len(intervals) >= CAP:
            break

    return intervals

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()