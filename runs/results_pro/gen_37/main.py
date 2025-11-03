# EVOLVE-BLOCK-START

def construct_intervals(depth=None):
  """
  Enhanced cycle of seeds, bridges, and scaling to push FirstFit colors higher
  while keeping omega moderate. This version expands the pattern banks
  and introduces an additional seed and longer scaling schedule.
  """

  if depth is None:
    depth = 4

  # Expanded pattern banks (8 patterns each) for richer diversification
  START_SETS = [
    [2, 6, 10, 14], [3, 7, 11, 15], [4, 8, 12, 16], [5, 9, 13, 17],
    [6, 10, 14, 18], [7, 11, 15, 19], [8, 12, 16, 20], [9, 13, 17, 21]
  ]
  BRIDGE_SETS = [
    [(1, 5), (12, 16), (4, 9), (8, 13)],
    [(2, 6), (13, 17), (5, 10), (9, 14)],
    [(3, 7), (14, 18), (6, 11), (10, 15)],
    [(4, 8), (15, 19), (7, 12), (11, 16)],
    [(5, 9), (16, 20), (8, 13), (12, 17)],
    [(6, 10), (17, 21), (9, 14), (13, 18)],
    [(7, 11), (18, 22), (10, 15), (14, 19)],
    [(8, 12), (19, 23), (11, 16), (15, 20)],
  ]

  # Slightly longer and more varied scaling factors
  SCALING_LIST = [1.0, 1.25, 0.88, 1.18, 1.05, 0.92, 1.12, 0.85]

  # Start with a two-spine seed and a small central spine for greater diversity
  T = [(0, 1), (2, 3), (6, 7)]

  for k in range(depth):
    lo = min(l for l, r in T)
    hi = max(r for l, r in T)
    gamma = SCALING_LIST[k % len(SCALING_LIST)]
    delta = (hi - lo) * gamma
    S = []
    start_list = START_SETS[k % len(START_SETS)]
    for start in start_list:
      S += [(delta * start + l - lo, delta * start + r - lo) for l, r in T]
    bridge_set = BRIDGE_SETS[k % len(BRIDGE_SETS)]
    S += [ (delta * a, delta * b) for (a, b) in bridge_set ]
    T = S
  return T

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()