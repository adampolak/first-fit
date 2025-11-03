# EVOLVE-BLOCK-START

def construct_intervals():
  """
  Construct a sequence of intervals of real line,
  in the order in which they are presented to FirstFit,
  so that it maximizes the number of colors used by FirstFit
  divided by the maximum number of intervals that cover a single point

  This version deepens the recursive cap construction to four layers and
  alternates start patterns across layers to increase asymmetry:
  (2,6,10,14) on even layers and (3,7,11,15) on odd layers.
  The canonical bridge intervals from Figure 4 are kept each layer.

  Returns:
    intervals: list of tuples, each tuple (l, r) represents an open interval from l to r
  """

  T = [(0, 1)]
  # Four recursive expansions with alternating start patterns
  for i in range(4):
    lo = min(l for l, r in T)
    hi = max(r for l, r in T)
    delta = hi - lo
    S = []
    starts = (2, 6, 10, 14) if (i % 2 == 0) else (3, 7, 11, 15)
    for start in starts:
      S += [(delta * start + l - lo, delta * start + r - lo) for l, r in T]
    # Canonical bridges
    S += [
      (delta * 1, delta * 5),
      (delta * 12, delta * 16),
      (delta * 4, delta * 9),
      (delta * 8, delta * 13),
    ]
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