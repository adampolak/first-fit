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

  # Use the classic fixed start-pattern to preserve strong KT coupling
  starts = (2, 6, 10, 14)
  T = [(0, 1)]
  for round_idx in range(6):
    lo = min(l for l, r in T)
    hi = max(r for l, r in T)
    delta = hi - lo
    # build translated blocks in sequence (no interleaving)
    blocks = []
    for start in starts:
      blocks.append([(delta * start + l - lo, delta * start + r - lo) for l, r in T])
    S = []
    for blk in blocks:
      S.extend(blk)
    # connectors based on the fixed starts
    s0, s1, s2, s3 = starts
    connectors = [
      (delta * (s0 - 1), delta * (s1 - 1)),  # left cap
      (delta * (s2 + 2), delta * (s3 + 2)),  # right cap
      (delta * (s0 + 2), delta * (s2 - 1)),  # cross 1
      (delta * (s1 + 2), delta * (s3 - 1)),  # cross 2
    ]
    for a, b in connectors:
      S.append((a, b))
    T = S
  # No micro-phase: returning the classic KT construction maximizes FF pressure here
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