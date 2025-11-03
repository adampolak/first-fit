# EVOLVE-BLOCK-START

def construct_intervals(rounds=3, seed_lo=0.0, seed_scale=1.0, extra_rounds=0):
    """
    Hierarchical shrink + spine-blocker construction.

    Changes summary:
    - Normalize the pattern each round and apply an explicit shrinking factor to control interval widths.
    - Place scaled copies at widely separated offsets (block_spacing) so global omega stays small.
    - Insert a small number of long "spine" blockers first in each block (to occupy small colors),
      then present short normalized intervals; this arrival-order increases FirstFit pressure.
    - Rotate among a deterministic bank of 4-start patterns and templates; add small deterministic jitter.
    - Optional extra_rounds reuse the same shrinking principle with denser starts.
    """
    # Seed with several disjoint unit intervals to promote early overlap coupling
    T = [
        (seed_lo, seed_lo + 1.0 * seed_scale),
        (seed_lo + 2.0 * seed_scale, seed_lo + 3.0 * seed_scale),
        (seed_lo + 4.0 * seed_scale, seed_lo + 5.0 * seed_scale),
        (seed_lo + 6.0 * seed_scale, seed_lo + 7.0 * seed_scale),
    ]

    # A small bank of four-interval gadgets (templates) — used at reduced scale each round
    templates = [
        [(1, 5), (12, 16), (4, 9), (8, 13)],
        [(0.5, 4.5), (11, 15), (3.5, 8.5), (7, 12)],
        [(1, 4), (6, 9), (3, 7), (9, 13)],
        [(2, 6), (7, 11), (0, 3), (10, 14)]
    ]

    # Deterministic rotation of compact 4-start patterns to vary alignment across rounds
    start_patterns = [
        (2, 6, 10, 14),
        (1, 5, 9, 13),
        (3, 7, 11, 15),
        (2, 4, 8, 12),
    ]

    base_block_spacing = 20.0   # base separation between blocks (keeps blocks mostly disjoint)
    base_shrink = 0.5           # per-round geometric shrink (width_scale = base_shrink ** round_idx)
    spine_count = 3             # number of long blockers per block (keeps omega small)

    for round_idx in range(rounds):
        lo = min(l for l, r in T)
        hi = max(r for l, r in T)
        delta = hi - lo
        if delta <= 0:
            delta = 1.0

        # normalize existing pattern into [0,1] before scaling/placing
        normT = [((l - lo) / delta, (r - lo) / delta) for (l, r) in T]

        # pick deterministic pattern/template for this round
        starts = list(start_patterns[round_idx % len(start_patterns)])
        template = templates[round_idx % len(templates)]

        # width scale shrinks with depth (bounded below to avoid numerical underflow)
        width_scale = max(0.12, (base_shrink ** round_idx))
        # mild growth of spacing to ensure distant rounds don't accidentally interact
        block_spacing = base_block_spacing * (1.1 ** round_idx)

        S = []

        # Build blocks; within each block add long spine blockers first (to occupy small colors),
        # then place the short, scaled normalized intervals.
        for idx, start in enumerate(starts):
            # deterministic tiny jitter to break symmetry
            jitter = ((round_idx * 7 + idx * 13) % 17) * 0.01
            offset = start * block_spacing + jitter

            # Add spine blockers first (occupy the small colors)
            for s in range(spine_count):
                s_lo = offset - 0.25 + s * 0.01
                s_hi = offset + width_scale + 0.25
                S.append((s_lo, s_hi))

            # Add scaled normalized intervals (alternate inner order to reduce color reuse)
            base_norm = normT if (idx % 2 == 0) else list(reversed(normT))
            for (nl, nr) in base_norm:
                S.append((offset + width_scale * nl, offset + width_scale * nr))

        # light bridges between adjacent blocks to couple colors (but avoid creating huge cliques)
        for i in range(len(starts) - 1):
            off_i = starts[i] * block_spacing
            off_j = starts[i + 1] * block_spacing
            mid1 = off_i + 0.5 * width_scale
            mid2 = off_j + 0.5 * width_scale
            S.append((mid1 + 0.1, mid2 - 0.1))

        # place the chosen template gadgets scaled down and located near the central region
        central_offset = sum(starts) * block_spacing / (len(starts) * 1.0)
        for (a, b) in template:
            # scale template coordinates to fit inside the central block region
            t_lo = central_offset + width_scale * (a / 4.0)
            t_hi = central_offset + width_scale * (b / 4.0)
            S.append((t_lo, t_hi))

        # short caps near each block center to pressure FirstFit locally
        for start in starts:
            off = start * block_spacing
            for k in range(2):
                cap_mid = off + (0.2 + 0.1 * k) * width_scale
                S.append((cap_mid - 0.05 * width_scale, cap_mid + 0.05 * width_scale))

        T = S

    # Optional extra shrinking rounds with a denser set of starts (still bounded width_scale)
    for er in range(extra_rounds):
        lo = min(l for l, r in T)
        hi = max(r for l, r in T)
        delta = hi - lo
        if delta <= 0:
            delta = 1.0
        normT = [((l - lo) / delta, (r - lo) / delta) for (l, r) in T]

        extended_starts = [2, 4, 6, 8, 12, 16, 20, 24]
        width_scale = max(0.08, (base_shrink ** (rounds + er)))
        block_spacing = base_block_spacing * (1.1 ** (rounds + er))

        S = []
        for idx, start in enumerate(extended_starts):
            jitter = ((er * 11 + idx * 19) % 23) * 0.01
            offset = start * block_spacing + jitter
            # spine first
            for s in range(spine_count):
                S.append((offset - 0.25 + s * 0.01, offset + width_scale + 0.25))
            base_norm = normT if (idx % 2 == 0) else list(reversed(normT))
            for (nl, nr) in base_norm:
                S.append((offset + width_scale * nl, offset + width_scale * nr))

        # add a sparse template in the center of the extended pattern
        center_idx = len(extended_starts) // 2
        center_offset = extended_starts[center_idx] * block_spacing
        for (a, b) in templates[0]:
            S.append((center_offset + width_scale * (a / 4.0),
                      center_offset + width_scale * (b / 4.0)))
        T = S

    return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()