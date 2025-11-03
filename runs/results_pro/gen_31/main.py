# EVOLVE-BLOCK-START

def construct_intervals():
    """
    Hybrid gadget: four-layer recursive expansion with alternating shift patterns,
    plus canonical bridge intervals. Outputs a sequence of open intervals.

    Improvements applied:
    - cycle a set of four start-factor patterns to break symmetry;
    - cycle bridge configurations across layers;
    - apply a small gamma-scaling schedule to the span each layer to diversify geometry;
    - use two disjoint seed intervals to break early-stage symmetry;
    - interleave copies when expanding a layer so FirstFit sees a more entangled order.
    """

    def expand_layer(template, starts, bridges, gamma=1.0, interleave=True, layer_index=0):
        """
        Given a list of intervals `template`, produce one expansion layer:
        - replicate `template` shifted by each factor in `starts` using an effective span
          (span * gamma)
        - optionally interleave the replicated copies to increase color mixing
        - add bridge intervals scaled by the effective span (split before/after interleaving)
        """
        lo = min(l for l, r in template)
        hi = max(r for l, r in template)
        span = hi - lo
        eff = span * gamma if span > 0 else 1.0
        # build shifted copies
        copies = []
        for s in starts:
            offset = eff * s - lo
            shifted = [(l + offset, r + offset) for (l, r) in template]
            copies.append(shifted)
        result = []
        if interleave and template:
            half = len(bridges) // 2
            # put some bridges before to break regularity on even layers
            if (layer_index % 2) == 0:
                for a, b in bridges[:half]:
                    result.append((eff * a, eff * b))
            # interleave entries by original template index
            n = len(template)
            for k in range(n):
                for copy in copies:
                    result.append(copy[k])
            # remaining bridges
            for a, b in bridges[half:]:
                result.append((eff * a, eff * b))
        else:
            # blockwise append copies then all bridges
            for copy in copies:
                result.extend(copy)
            for a, b in bridges:
                result.append((eff * a, eff * b))
        return result

    def build(depth):
        """
        Build the gadget by repeatedly expanding the template `depth` times.
        We cycle several start-factor patterns, bridge variants and a small gamma schedule.
        """
        # initial seed: use two disjoint unit intervals to break symmetry
        T = [(0.0, 1.0), (3.0, 4.0)]
        # four start-factor patterns cycled across layers
        patterns = [
            (2, 6, 10, 14),
            (3, 7, 11, 15),
            (4, 8, 12, 16),
            (5, 9, 13, 17),
        ]
        # bridge configurations (rotations/permutations of the canonical quartet)
        bridge_variants = [
            [(1,5), (12,16), (4,9), (8,13)],
            [(1,5), (4,9), (12,16), (8,13)],
            [(4,9), (1,5), (8,13), (12,16)],
            [(8,13), (4,9), (12,16), (1,5)],
        ]
        # gentle gamma schedule to vary effective spacing
        gamma_seq = [1.00, 1.05, 0.98, 1.10]
        for i in range(depth):
            starts = patterns[i % len(patterns)]
            bridges = bridge_variants[i % len(bridge_variants)]
            gamma = gamma_seq[i % len(gamma_seq)]
            T = expand_layer(T, starts, bridges, gamma=gamma, interleave=True, layer_index=i)
        return T

    def normalize(intervals, scale=10):
        """
        Shift all intervals to positive, multiply by `scale`, and
        round to integers so that open intervals remain valid.
        """
        min_l = min(l for l, r in intervals)
        shift = (-min_l + 1.0) if min_l <= 0.0 else 0.0
        normalized = []
        for (l, r) in intervals:
            L = int((l + shift) * scale)
            R = int((r + shift) * scale + 0.9999)
            if R <= L:
                R = L + 1
            normalized.append((float(L), float(R)))
        return normalized

    # Build with 4 recursive layers to push FF >> opt (diversified)
    raw = build(depth=4)
    # Map to integer grid and clean up
    return normalize(raw)

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()