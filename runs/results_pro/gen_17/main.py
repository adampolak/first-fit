# EVOLVE-BLOCK-START

def construct_intervals():
    """
    Hybrid gadget: four-layer recursive expansion with alternating shift patterns,
    plus canonical bridge intervals. Outputs a sequence of open intervals.
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
            for (l, r) in template:
                result.append((l + offset, r + offset))
        # add bridges (each specified as (a,b) factors)
        for a, b in bridges:
            result.append((span * a, span * b))
        return result

    def build(depth):
        """
        Build the gadget by repeatedly expanding the template `depth` times.
        We alternate two start‐factor patterns to induce asymmetry.
        """
        # initial seed: single unit interval
        T = [(0.0, 1.0)]
        # fixed bridge-interval factors from Figure 4
        bridge_factors = [(1,5), (12,16), (4,9), (8,13)]
        # two alternating start-factor patterns
        patterns = [
            (2,6,10,14),
            (3,7,11,15),
        ]
        for i in range(depth):
            starts = patterns[i % 2]
            T = expand_layer(T, starts, bridge_factors)
        return T

    def normalize(intervals, scale=15):
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

    # Build with 5 recursive layers to push FF further while keeping omega modest
    raw = build(depth=5)
    # Map to integer grid and clean up
    # Also append a secondary gadget translated by a fixed gap to increase interaction
    min_l = min(l for l, r in raw)
    max_r = max(r for l, r in raw)
    span = max_r - min_l if max_r > min_l else 1.0
    gap = span * 0.4
    second = [(l + gap, r + gap) for (l, r) in raw]
    combined = raw + second
    return normalize(combined, scale=15)

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()