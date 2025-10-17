# EVOLVE-BLOCK-START

def construct_intervals():
  """
  Construct a sequence of intervals of real line,
  in the order in which they are presented to FirstFit,
  so that it maximizes the number of colors used by FirstFit
  divided by the maximum number of intervals that cover a single point
  
  The initial implementation uses the construction from
  Figure 2 in https://arxiv.org/abs/1506.00192
  
  Returns:
    intervals: list of tuples, each tuple (l, r) represents an open interval from l to r
  """
  return [(1, 4), (7, 10), (3, 6), (5, 8)]

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()