# EVOLVE-BLOCK-START

def construct_intervals():
  """
  Construct a sequence of intervals of real line,
  in the order in which they are presented to FirstFit,
  so that it tends to maximize the number of colors used by FirstFit
  divided by the maximum number of intervals that cover a single point.

  New approach: rotating_template_wave
  - Depth-controlled multi-block expansion.
  - A bank of templates (offset patterns) rotated deterministically per level.
  - Each level places multiple translated copies of the current gadget with alternating orientation.
  - Adds local connectors between blocks to propagate color usage but keeps omega small.
  - Deterministic and keeps total number of intervals under practical limits.

  Returns:
    intervals: list of tuples, each tuple (l, r) represents an open interval from l to r
  """

  # Depth controls how many recursive expansion rounds we perform.
  depth = 4  # tuned to stay well under 10000 intervals
  # Deterministic bank of templates. Each template is a list of starting offsets
  # relative to the current span. We apply them in order, flipping orientation for every other block.
  templates = [
    [2, 6, 10, 14, 18],  # template A
    [3, 7, 11, 15, 19],  # template B
    [1, 5, 9, 13, 17],   # template C
    [4, 8, 12, 16, 20],  # template D
    [2, 5, 9, 14, 19],   # template E
  ]

  # Start with a single base gadget
  T = [(0, 1)]
  for depth_idx in range(depth):
    lo = min(l for l, r in T)
    hi = max(r for l, r in T)
    span = hi - lo
    S = []

    # Choose a deterministic template for this level
    template = templates[depth_idx % len(templates)]
    # Place copies for each start in the chosen template
    for block_idx, start in enumerate(template):
      # Alternate orientation of the base gadget to break symmetry
      block_T = T[::-1] if (block_idx % 2 == 1) else T
      for (l, r) in block_T:
        S.append((span * start + l - lo, span * start + r - lo))

    # Connect adjacent blocks to propagate colors across blocks
    # These connectors are kept modest to avoid exploding omega
    for i in range(len(template) - 1):
      a = span * (1 + 4 * i)
      b = span * (5 + 4 * i)
      S.append((a, b))

    # Cross-connectors to couple non-adjacent blocks without increasing clique
    for i in range(len(template) - 2):
      a = span * (4 + 4 * i)
      b = span * (9 + 4 * i)
      S.append((a, b))

    # Small additional local intervals to seed color pressure
    S.append((span * 2, span * 3))
    S.append((span * 7, span * 8))

    T = S

  return T

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()
# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()