# EVOLVE-BLOCK-START

def construct_intervals():
  """
  Construct a sequence of intervals of real line,
  in the order in which they are presented to FirstFit,
  so that it maximizes the number of colors used by FirstFit
  divided by the maximum number of intervals that cover a single point

  This implementation follows the recursive "cap" construction
  inspired by Figure 4 in https://arxiv.org/abs/1506.00192 and
  increases the recursion depth (subject to a size budget) to amplify
  the FF/omega ratio. After d expansions, FF uses ≈ 3d+1 colors
  while omega is d+1, so the ratio approaches 3 from below.

  Returns:
    intervals: list of tuples, each tuple (l, r) represents an open interval from l to r
  """

  # Choose the deepest recursion allowed by an interval-count budget
  # to improve FF/omega while keeping the instance compact.
  def _count_intervals(d: int) -> int:
    # Recurrence: n_d = 4 n_{d-1} + 4, n_0 = 1  =>  n_d = (7*4^d - 4)/3
    return (7 * (4 ** d) - 4) // 3

  max_intervals = 10_000
  target_depth = 6  # depth 6 gives ~9556 intervals and ratio ~19/7 ≈ 2.714
  depth = target_depth
  while depth > 1 and _count_intervals(depth) > max_intervals:
    depth -= 1

  # Cap construction
  T = [(0, 1)]
  for _ in range(depth):
    lo = min(l for l, r in T)
    hi = max(r for l, r in T)
    delta = hi - lo
    S = []
    for start in (2, 6, 10, 14):
      S += [(delta * start + l - lo, delta * start + r - lo) for l, r in T]
    S += [
      (delta * 1, delta * 5),
      (delta * 12, delta * 16),
      (delta * 4, delta * 9),
      (delta * 8, delta * 13)
    ]
    T = S

  # Normalize by shifting so the minimum left endpoint is 0.
  # This shrink step preserves all intersections for open intervals.
  base = min(l for l, r in T)
  if base != 0:
    T = [(l - base, r - base) for l, r in T]

  return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()