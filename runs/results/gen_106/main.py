# EVOLVE-BLOCK-START

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

  T = [(0, 1)]
  for i in range(4):
    lo = min(l for l, r in T)
    hi = max(r for l, r in T)
    delta = hi - lo
    # define blockers
    blockers = [
      (delta * 1, delta * 5),
      (delta * 12, delta * 16),
      (delta * 4, delta * 9),
      (delta * 8, delta * 13)
    ]
    # generate copies
    copies = []
    for start in (2, 6, 10, 14):
      copies += [(delta * start + l - lo, delta * start + r - lo) for l, r in T]
    # interleave blockers and copies based on iteration parity
    if i % 2 == 0:
      S = blockers + copies
    else:
      S = copies + blockers
    T = S
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