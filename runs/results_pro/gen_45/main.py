# EVOLVE-BLOCK-START

def construct_intervals():
    """
    Construct an interval sequence designed to maximize the FirstFit/omega ratio.
    Strategy:
    - Build a 4-layer Figure-4-inspired recursive gadget (cap construction) using
      alternating shift patterns and canonical bridge intervals. This backbone is
      known to produce alg≈13 with omega≈5 (ratio≈2.60).
    - Normalize to an integer grid and adversarially order by increasing length,
      a safe reordering that preserves the intersection graph but tends to
      increase FirstFit usage.

    The implementation borrows the modular layer-expansion style from the
    hybrid code, but uses conservative parameters (no gamma drift, no interleaving)
    to keep the strong 2.60 ratio backbone intact.
    Returns:
      List of (l, r) open-interval tuples with float endpoints.
    """

    import math

    # Canonical bridge-interval factors (Figure 4 gadget)
    BRIDGES = [(1, 5), (12, 16), (4, 9), (8, 13)]
    # Alternating start-factor patterns
    PATTERNS = [
        (2, 6, 10, 14),
        (3, 7, 11, 15),
    ]

    def expand_layer(template, starts, bridges, gamma=1.0, interleave=False, layer_index=0):
        """
        Expand one layer:
        - Replicate `template` shifted by each factor in `starts` using effective span eff = span*gamma.
        - Append `bridges` scaled by eff.
        - `interleave` is supported (from the hybrid code) but disabled by default to preserve
          the backbone behavior yielding the best ratio.
        """
        lo = min(l for l, r in template)
        hi = max(r for l, r in template)
        span = hi - lo
        eff = span * gamma if span > 0 else 1.0

        # Build shifted copies of the template
        copies = []
        for s in starts:
            offset = eff * s - lo
            shifted = [(l + offset, r + offset) for (l, r) in template]
            copies.append(shifted)

        result = []
        if interleave and template:
            # Interleave entries by original template index
            # Note: kept for completeness but not used (interleave=False).
            half = len(bridges) // 2
            if (layer_index % 2) == 0:
                for a, b in bridges[:half]:
                    result.append((eff * a, eff * b))
            n = len(template)
            for k in range(n):
                for copy in copies:
                    result.append(copy[k])
            for a, b in bridges[half:]:
                result.append((eff * a, eff * b))
        else:
            # Append blockwise copies then bridges (matches best-performing backbone behavior)
            for copy in copies:
                result.extend(copy)
            for a, b in bridges:
                result.append((eff * a, eff * b))
        return result

    def build_backbone(depth=4):
        """
        Build the cap gadget backbone across `depth` layers using alternating
        start patterns and canonical bridges.
        """
        T = [(0.0, 1.0)]  # seed
        for i in range(depth):
            starts = PATTERNS[i % 2]
            T = expand_layer(T, starts, BRIDGES, gamma=1.0, interleave=False, layer_index=i)
        return T

    def normalize_and_order(intervals, scale=10):
        """
        Shift to positive, scale to integer grid, and sort by increasing length (then left endpoint).
        Sorting changes only arrival order (not the intersection graph), often increasing FirstFit usage.
        """
        min_l = min(l for l, r in intervals)
        shift = (-min_l + 1.0) if min_l <= 0.0 else 0.0
        norm = []
        for (l, r) in intervals:
            L = int(math.floor((l + shift) * scale))
            R = int(math.ceil((r + shift) * scale))
            if R <= L:
                R = L + 1
            norm.append((float(L), float(R)))
        # Adversarial arrival order: shortest first
        norm.sort(key=lambda iv: (iv[1] - iv[0], iv[0], iv[1]))
        return norm

    # Build the proven backbone
    raw = build_backbone(depth=4)

    # Normalize and adversarially order
    final = normalize_and_order(raw, scale=10)

    # Deduplicate exactly-equal intervals while preserving order
    seen = set()
    dedup = []
    for iv in final:
        if iv in seen:
            continue
        seen.add(iv)
        dedup.append(iv)

    return dedup

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()