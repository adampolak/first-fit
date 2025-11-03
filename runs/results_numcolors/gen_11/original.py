# EVOLVE-BLOCK-START

def construct_intervals(rounds=2):
  """
  Construct a sequence of intervals of real line,
  in the order in which they are presented to FirstFit,
  so that it tends to maximize the number of colors used by FirstFit
  divided by the maximum number of intervals that cover a single point
  through a recursive wave construction.

  The construction generalizes the two-round pattern used in Figure 4 of
  https://arxiv.org/abs/1506.00192 to a multi-round recursive expansion.

  Returns:
    intervals: list of tuples, each tuple (l, r) represents an open interval from l to r
  """
  T = [(0, 1)]
  for round_idx in range(rounds):
    lo = min(l for l, r in T)
    hi = max(r for l, r in T)
    delta = hi - lo
    S = []
    # In earlier rounds we apply a compact set of shifts; in later rounds,
    # we widen the set of shifts to inject more waves and increase color usage.
    # Density-increasing starts: first round uses a core set, subsequent rounds add more starts
    if round_idx == 0:
      starts = [2, 4, 6, 8, 10, 12, 14, 16]
    else:
      starts = [2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 24, 28]
    for start in starts:
      S += [(delta * start + l - lo, delta * start + r - lo) for l, r in T]
    # additional local blocks to create overlapping waves without increasing omega too much
    S += [
      (delta * 1, delta * 5),
      (delta * 12, delta * 16),
      (delta * 4, delta * 9),
      (delta * 8, delta * 13),
      # extra caps to push color usage further without bloating omega
      (delta * 0.5, delta * 2.0),
      (delta * 14.0, delta * 15.0)
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