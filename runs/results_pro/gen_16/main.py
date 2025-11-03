# EVOLVE-BLOCK-START

def construct_intervals():
  """
  Construct a sequence of intervals on the real line presented to FirstFit,
  based on the Figure 4 gadget from https://arxiv.org/abs/1506.00192,
  iterated four times to further inflate the color usage while keeping the
  clique number near 3–4.

  Improvements:
  - Cycle through four start-factor patterns (rather than a 2-alternate),
    and cycle through a small family of bridge-sets. This increases asymmetry
    and cross-coupling between layers, a heuristic that tends to push FirstFit
    to open more colors without increasing the clique too much.
  - Normalize the final intervals to a compact integer grid.
  Returns:
    intervals: list of tuples (l, r) representing open intervals.
  """

  # Start with a single unit interval.
  T = [(0.0, 1.0)]

  # define a small cycle of start-factor patterns to vary placements across layers
  patterns = [
    (2, 6, 10, 14),
    (3, 7, 11, 15),
    (4, 8, 12, 16),
    (5, 9, 13, 17),
  ]

  # define a small cycle of bridge-sets to slightly change coupling each layer
  bridge_sets = [
    [(1, 5), (12, 16), (4, 9), (8, 13)],  # baseline
    [(2, 6), (11, 15), (5, 10), (9, 14)],  # shifted variant
    [(3, 7), (10, 14), (6, 11), (7, 12)],  # midband couplers
    [(4, 8), (13, 17), (7, 12), (8, 13)],  # heavier right cap
  ]

  # Perform four recursive expansions, cycling patterns and bridge sets
  depth = 4
  for i in range(depth):
    lo = min(l for l, r in T)
    hi = max(r for l, r in T)
    delta = hi - lo
    S = []
    starts = patterns[i % len(patterns)]
    bridges = bridge_sets[i % len(bridge_sets)]
    for start in starts:
      offset = delta * start - lo
      S += [(l + offset, r + offset) for (l, r) in T]
    # add the chosen bridge-set scaled by current span
    for a, b in bridges:
      S.append((delta * a, delta * b))
    T = S

  # Normalize into a compact integer grid to keep coordinates moderate for the evaluator.
  # Preserve ordering and open-interval semantics by mapping via linear scaling + shift.
  min_l = min(l for l, r in T)
  if min_l <= 0:
    shift = -min_l + 1.0
    T = [(l + shift, r + shift) for (l, r) in T]
  lo = min(l for l, r in T)
  hi = max(r for l, r in T)
  span = hi - lo if hi > lo else 1.0
  # target a modest span (about 2000) to avoid excessively large coordinates
  scale = max(1.0, 2000.0 / span)
  normalized = []
  for (l, r) in T:
    L = int((l - lo) * scale)
    R = int((r - lo) * scale + 0.9999)
    if R <= L:
      R = L + 1
    normalized.append((float(L), float(R)))

  return normalized

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()