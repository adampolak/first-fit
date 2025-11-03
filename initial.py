# EVOLVE-BLOCK-START

def construct3vs2:
  # The construction from Figure 2 in https://arxiv.org/abs/1506.00192
  return [(1, 4), (7, 10), (3, 6), (5, 8)]

def construct4vs2():
  # The construction from Figure 3 in https://arxiv.org/abs/1506.00192
  return [
    (2,3),
    (6,7),
    (10,11),
    (14,15),
    (1,5),
    (12,16),
    (4,9),
    (8,13),
  ]

def construct7vs4():
  # The construction from Figure 4 in https://arxiv.org/abs/1506.00192
  T = [(0, 1)]
  for _ in range(2):
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
  return T

def construct_intervals():
  """
  Construct a sequence of intervals of real line,
  in the order in which they are presented to FirstFit,
  so that it maximizes the number of colors used by FirstFit
  divided by the maximum number of intervals that cover a single point
  
  The initial implementation uses the construction from
  Figure 4 in https://arxiv.org/abs/1506.00192
  
  Returns:
    intervals: list of tuples, each tuple (l, r) represents an open interval from l to r
  """
  return construct7vs4()
  
# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()