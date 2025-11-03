# EVOLVE-BLOCK-START

def construct_intervals(iterations=4):
  """
  Construct a sequence of intervals of real line,
  in the order in which they are presented to FirstFit,
  so that it maximizes the number of colors used by FirstFit
  divided by the maximum number of intervals that cover a single point

  The implementation uses the construction from
  Figure 4 in https://arxiv.org/abs/1506.00192 and increases
  the number of recursive expansions to produce a stronger adversary.

  Arguments:
    iterations: number of recursive expansion steps (default 4)

  Returns:
    intervals: list of tuples, each tuple (l, r) represents an open interval from l to r
  """

  T = [(0, 1)]
  # Each iteration replaces T by four translated/scaled copies plus four
  # longer intervals that connect the copies in the manner of Figure 4.
  for i in range(iterations):
    lo = min(l for l, r in T)
    hi = max(r for l, r in T)
    delta = hi - lo
    S = []
    # dynamic offsets: extra copy on first iteration to diversify branching
    if i == 0:
      offsets = (2, 6, 10, 14, 18)
    else:
      offsets = (2, 6, 10, 14)
    for start in offsets:
      S += [(delta * start + l - lo, delta * start + r - lo) for l, r in T]
    # Prepare the four connecting long intervals (blockers)
    blockers = [
      (delta * 1,  delta * 5),
      (delta * 12, delta * 16),
      (delta * 4,  delta * 9),
      (delta * 8,  delta * 13)
    ]
    # Alternate insertion schedule: even iterations after copies, odd before
    if i % 2 == 0:
      # even: append blockers
      S += blockers
    else:
      # odd: prepend blockers
      S = blockers + S
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