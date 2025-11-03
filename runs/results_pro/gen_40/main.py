# EVOLVE-BLOCK-START

def construct_intervals(depth=None):
  """
  Construct a sequence of intervals of real line,
  in the order in which they are presented to FirstFit,
  so that it maximizes the number of colors used by FirstFit
  divided by the maximum number of intervals that cover a single point

  This implementation follows the recursive "cap" construction
  inspired by Figure 4 in https://arxiv.org/abs/1506.00192 and
  increases the recursion depth to amplify the FF/omega ratio.

  Returns:
    intervals: list of tuples, each tuple (l, r) represents an open interval from l to r
  """

  # Allow explicit depth control to explore smarter configurations.
  # If depth is not provided, fall back to the traditional depth=4.
  if depth is None:
    depth = 5

  # Pattern banks to introduce variety across recursion levels.
  START_SETS = [
    [2, 6, 10, 14],
    [3, 7, 11, 15],
    [4, 8, 12, 16],
    [5, 9, 13, 17],
  ]
  BRIDGE_SETS = [
    [(1, 5), (12, 16), (4, 9), (8, 13)],
    [(2, 6), (13, 17), (5, 10), (9, 14)],
    [(3, 7), (14, 18), (6, 11), (10, 15)],
    [(4, 8), (15, 19), (7, 12), (11, 16)],
  ]

  SCALING_LIST = [1.0, 1.2, 0.6, 1.15, 1.05]
  # initial seeds: three disjoint spine intervals to diversify recursive growth
  T = [(0, 1), (2, 3), (4, 5)]
  for k in range(depth):
    lo = min(l for l, r in T)
    hi = max(r for l, r in T)
    gamma = SCALING_LIST[k % len(SCALING_LIST)]
    delta = (hi - lo) * gamma
    S = []
    # Use a rotating start set to diversify the expanded copies of T
    start_list = START_SETS[k % len(START_SETS)]
    for start in start_list:
      S += [(delta * start + l - lo, delta * start + r - lo) for l, r in T]
    # Use a rotating bridge pattern to insert cross-links without inflating omega
    bridge_set = BRIDGE_SETS[k % len(BRIDGE_SETS)]
    S += [ (delta * a, delta * b) for (a, b) in bridge_set ]
    T = S
  # Add two global caps to couple branches across the gadget
  lo = min(l for l, r in T)
  hi = max(r for l, r in T)
  span = hi - lo
  # caps spanning first and last translations of the first two start patterns
  T.append((lo + span * START_SETS[0][0], lo + span * START_SETS[0][-1]))
  T.append((lo + span * START_SETS[1][0], lo + span * START_SETS[1][-1]))
  return T

  # return [  # Figure 3, OPT=2, FF=4
  #   (2,3),
  #   (6,7),
  #   (10,11),
  #   (14,15),
  #   (1,5),
  #   (12,16),
  #   (4,9),
  #   (8,13),
  # ]

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()