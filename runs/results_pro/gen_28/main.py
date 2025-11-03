# EVOLVE-BLOCK-START

def construct_intervals():
    """
    Construct a sequence of open intervals presented to FirstFit that aims to
    maximize the ratio (#colors used by FirstFit) / (clique number = offline optimum).
    Structural redesign:
      - Four-start translation pattern cycles across iterations (A,B,C,D).
      - Bridge patterns also cycle across iterations (A,B,C,D).
      - Deterministic gamma schedule per layer diversifies span growth.
      - Two-seed start to break early symmetry.
    Returns:
      list[(float,float)]: Open intervals in arrival order.
    """
    import math

    # ---------- Pattern cycles (deterministic) ----------
    # Start-shift cycles (four 4-tuples), scaled each layer by span*gamma_i
    start_cycles = [
        (2, 6, 10, 14),  # A
        (3, 7, 11, 15),  # B
        (4, 8, 12, 16),  # C
        (5, 9, 13, 17),  # D
    ]

    # Bridge cycles (four 4-tuples), all scaled by span*gamma_i each layer
    # These sets are chosen to couple far and near blocks in different pairings across layers.
    bridge_cycles = [
        # A
        [(1, 5), (12, 16), (4, 9), (8, 13)],
        # B
        [(2, 6), (11, 15), (5, 10), (9, 14)],
        # C
        [(3, 7), (10, 14), (6, 11), (7, 12)],
        # D
        [(4, 8), (13, 17), (7, 12), (8, 13)],
    ]

    # Deterministic gamma schedule for span diversification
    gamma_seq = [1.00, 1.25, 0.95, 1.15]

    # ---------- Core utilities ----------
    def span_of(intervals):
        lo = min(l for l, r in intervals)
        hi = max(r for l, r in intervals)
        return lo, hi, hi - lo

    def replicate(template, starts, bridges, gamma):
        """
        Replicate 'template' four times using starts and inject four bridges.
        Each layer is uniformly scaled by `gamma` to diversify spans deterministically.
        """
        lo, hi, span = span_of(template)
        scaled = []
        # Four translated, scaled copies of template
        for s in starts:
            offset = (span * s - lo) * gamma
            for (l, r) in template:
                scaled.append((l * gamma + offset, r * gamma + offset))
        # Four bridges scaled by (span * gamma) and anchored at 0
        for (a, b) in bridges:
            scaled.append((span * a * gamma, span * b * gamma))
        return scaled

    def normalize(intervals, target_scale=10.0):
        # Shift to positive and map to integer grid while preserving openness.
        lo = min(l for l, r in intervals)
        if lo <= 0.0:
            shift = -lo + 1.0
            intervals = [(l + shift, r + shift) for (l, r) in intervals]
        # Compactify if the span is huge/small
        lo, hi, sp = span_of(intervals)
        scale = target_scale if sp == 0 else max(1.0, target_scale)
        # Snap to integer grid (open intervals)
        out = []
        for (l, r) in intervals:
            L = int(math.floor(l * scale))
            R = int(math.ceil(r * scale))
            if R <= L:
                R = L + 1
            out.append((float(L), float(R)))
        # Deduplicate while preserving order
        seen = set()
        final = []
        for iv in out:
            if iv not in seen:
                seen.add(iv)
                final.append(iv)
        return final

    # ---------- Builder pipeline ----------
    def build_backbone(depth=4, use_two_seeds=True):
        """
        Build the multi-layer cap/backbone with cycling patterns and gamma schedule.
        Depth 4 remains a sweet spot; two seeds reduce early symmetry.
        """
        # Two disjoint seeds (optional) to diversify early FF color couplings.
        if use_two_seeds:
            T = [(0.0, 1.0), (2.0, 3.0)]
        else:
            T = [(0.0, 1.0)]

        for i in range(depth):
            starts = start_cycles[i % 4]
            bridges = bridge_cycles[i % 4]
            gamma = gamma_seq[i % 4]
            T = replicate(T, starts, bridges, gamma)
        return T

    # Build backbone with pattern/bridge/gamma cycling.
    backbone = build_backbone(depth=4, use_two_seeds=True)

    # Final normalization (compact mapping to integer grid)
    return normalize(backbone, target_scale=10.0)

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()