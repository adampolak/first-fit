# EVOLVE-BLOCK-START

def construct_intervals():
    """
    Construct an interval sequence that forces FirstFit to use a larger
    number of colors while keeping the clique number modest.
    This version uses a four-seed initialization, four alternating start-pattern
    layer expansions, and per-layer bridge scaling to inject diversity while
    preserving self-similarity and small omega.
    Returns:
      List of (l, r) open-interval tuples.
    """

    # Deterministic helpers and small utilities
    def expand_layer(template, starts, bridges, layer_idx, scale_seq):
        """
        Expand the current template once:
        - replicate the template translated by each factor in `starts`
        - add bridge intervals scaled by the current span and per-layer scale
        """
        lo = min(l for l, r in template)
        hi = max(r for l, r in template)
        span = hi - lo
        result = []
        # replicate by starts
        for s in starts:
            offset = span * s - lo
            for (l, r) in template:
                result.append((l + offset, r + offset))
        # bridges with per-layer scaling
        scale = scale_seq[layer_idx % len(scale_seq)]
        for a, b in bridges:
            result.append((span * a * scale, span * b * scale))
        return result

    def build(depth=4):
        """
        Build the gadget by repeatedly expanding the template `depth` times.
        We alternate four start-patterns to induce asymmetry across layers.
        """
        # initial four seeds (kept disjoint to keep omega small in early stages)
        T = [
            (0.0, 1.0),   # seed 0
            (5.0, 6.0),   # seed 1
            (11.0, 12.0), # seed 2
            (17.0, 18.0)  # seed 3
        ]
        # fixed bridge-interval factors from Figure 4 (canonical)
        bridge_factors = [(1,5), (12,16), (4,9), (8,13)]
        # four start-pattern cycles
        patterns = [
            (2, 6, 10, 14),
            (3, 7, 11, 15),
            (4, 8, 12, 16),
            (5, 9, 13, 17),
        ]
        scale_seq = [1.0, 1.25, 0.95, 1.15]  # deterministic per-layer scale variation

        for i in range(depth):
            starts = patterns[i % len(patterns)]
            T = expand_layer(T, starts, bridge_factors, i, scale_seq)
        return T

    def normalize(intervals, scale=10.0):
        """
        Normalize to positive coordinates and map to a compact integer grid.
        - shift so all l > 0
        - scale to target grid width
        - ensure open intervals remain valid
        """
        min_l = min(l for l, r in intervals)
        shift = (-min_l + 1.0) if min_l <= 0.0 else 0.0
        shifted = [(l + shift, r + shift) for (l, r) in intervals]
        lo = min(l for l, r in shifted)
        hi = max(r for l, r in shifted)
        span = hi - lo if hi > lo else 1.0
        # target around 2000 units across
        scale_eff = max(1.0, (2000.0 / span))
        normalized = []
        for (l, r) in shifted:
            L = int((l - lo) * scale_eff)
            R = int((r - lo) * scale_eff + 0.999999)
            if R <= L:
                R = L + 1
            normalized.append((float(L), float(R)))
        return normalized

    # Build backbone gadget with four layers
    raw = build(depth=4)

    # Normalize into a compact integer grid
    normalized = normalize(raw, scale=10.0)

    # Deduplicate exactly-equal intervals while preserving order
    seen = set()
    final = []
    for iv in normalized:
        if iv in seen:
            continue
        seen.add(iv)
        final.append(iv)

    return final

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()