# EVOLVE-BLOCK-START

def construct_intervals():
  """
  Construct a sequence of intervals on the real line presented to FirstFit,
  based on the Figure 4 gadget from https://arxiv.org/abs/1506.00192,
  iterated three times to force FirstFit to use many colors while
  the clique number (optimum) remains 3.
  Returns:
    intervals: list of tuples (l, r) representing open intervals.
  """

  # Start with a single unit interval.
  T = [(0.0, 1.0)]
  # Perform three recursive expansions (was two originally).
  for _ in range(3):
    # Compute current span
    lo = min(l for l, r in T)
    hi = max(r for l, r in T)
    delta = hi - lo
    S = []
    # Create four scaled copies of T at shifted positions
    for start in (2, 6, 10, 14):
      offset = delta * start - lo
      S += [ (l + offset, r + offset) for (l, r) in T ]
    # Add the four bridging intervals as in Figure 4
    S += [
      (delta * 1,  delta * 5),
      (delta * 12, delta * 16),
      (delta * 4,  delta * 9),
      (delta * 8,  delta * 13),
    ]
    T = S
  return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()