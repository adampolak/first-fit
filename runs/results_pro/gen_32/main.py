# EVOLVE-BLOCK-START

def construct_intervals():
    """
    Construct an interval sequence that forces FirstFit to use many colors
    while keeping the clique number small. We use a 5‐layer recursive
    Figure‑4 cap expansion with alternating shift patterns and canonical
    bridge intervals. Depth 5 empirically/analytically improves the
    FirstFit/omega ratio over depth 4 for this gadget class.
    Returns:
      List of (l, r) open‐interval tuples.
    """

    def expand_layer(template, starts, bridges):
        """
        Given a list of intervals `template`, produce one expansion layer:
        - replicate `template` shifted by each factor in `starts`
        - add each bridge interval scaled by the current span
        """
        lo = min(l for l, r in template)
        hi = max(r for l, r in template)
        span = hi - lo
        result = []
        # replicate
        for s in starts:
            offset = span * s - lo
            # append in the same order to preserve the FF forcing structure
            result.extend((l + offset, r + offset) for (l, r) in template)
        # add bridges (each specified as (a,b) factors)
        result.extend((span * a, span * b) for (a, b) in bridges)
        return result

    def build_multilayer(depth):
        """
        Build the gadget by repeatedly expanding the template `depth` times.
        We alternate two start‐shift patterns to induce asymmetry while
        preserving the classic cap invariants.
        """
        # initial seed: single unit interval
        T = [(0.0, 1.0)]
        # fixed bridge‐interval factors from Figure 4
        bridge_factors = [(1, 5), (12, 16), (4, 9), (8, 13)]
        # two alternating start‐factor patterns
        patterns = [
            (2, 6, 10, 14),
            (3, 7, 11, 15),
        ]
        for i in range(depth):
            starts = patterns[i % 2]
            T = expand_layer(T, starts, bridge_factors)
        return T

    def normalize(intervals, scale=10):
        """
        Shift all intervals to positive, multiply by `scale`, and
        round to integers so that open intervals remain valid.
        """
        # shift so all l > 0
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

    # Increase to 5 recursive layers to improve the FF/omega ratio.
    raw = build_multilayer(depth=5)
    # Map to integer grid and clean up
    return normalize(raw)

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()