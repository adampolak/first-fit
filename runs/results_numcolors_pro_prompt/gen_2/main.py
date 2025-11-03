# EVOLVE-BLOCK-START

def construct_intervals():
  """
  Construct a sequence of intervals of real line,
  in the order in which they are presented to FirstFit,
  so that it maximizes the number of colors used by FirstFit
  divided by the maximum number of intervals that cover a single point

  The initial implementation uses the construction from
  Figure 4 in https://arxiv.org/abs/1506.00192.  Here we
  dynamically increase the number of recursive replications
  to amplify the FirstFit color usage while keeping the total
  number of intervals under the practical limit.

  Returns:
    intervals: list of tuples, each tuple (l, r) represents an open interval from l to r
  """

  T = [(0, 1)]

  # Choose the maximal number of recursive replications so that
  # the total number of intervals remains under the allowed limit.
  max_intervals = 9999
  iter_count = 0
  curr_len = len(T)
  # Each iteration follows the recurrence next_len = 4*curr_len + 4
  while True:
    next_len = 4 * curr_len + 4
    if next_len > max_intervals:
      break
    curr_len = next_len
    iter_count += 1

  # Make sure we do at least the baseline number of iterations
  if iter_count < 2:
    iter_count = 2

  for _ in range(iter_count):
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