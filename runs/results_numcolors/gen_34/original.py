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

  # parameterize recursion depth and branching factor for stronger blow-up
  depth = 4
  branching = 7
  T = [(0, 1)]
  for _ in range(depth):
    lo = min(l for l, r in T)
    hi = max(r for l, r in T)
    delta = hi - lo
    S = []
    # insert scaled copies at offsets [2,6,10,...] based on branching
    starts = [2 + 4 * i for i in range(branching)]
    for idx, start in enumerate(starts):
      block_T = T[::-1] if (idx % 2 == 1) else T
      S += [(delta * start + l - lo, delta * start + r - lo) for l, r in block_T]
    # adjacent connectors to propagate FirstFit colors between neighbors
    for i in range(branching - 1):
      a = delta * (1 + 4 * i)
      b = delta * (5 + 4 * i)
      S.append((a, b))
    # cross connectors to couple non‐adjacent blocks without raising clique
    for i in range(branching - 2):
      a = delta * (4 + 4 * i)
      b = delta * (9 + 4 * i)
      S.append((a, b))
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