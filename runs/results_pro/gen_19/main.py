# EVOLVE-BLOCK-START

def construct_intervals(depth=None):
  """
  Construct a sequence of intervals of real line,
  in the order in which they are presented to FirstFit,
  so that it maximizes the number of colors used by FirstFit
  divided by the maximum number of intervals that cover a single point

  This implementation follows the recursive "cap" construction
  inspired by Figure 4 in https://arxiv.org/abs/1506.00192 and
  increases the recursion depth to amplify the FF/omega ratio.

  Returns:
    intervals: list of tuples, each tuple (l, r) represents an open interval from l to r
  """

  # Allow explicit depth control; explore slightly deeper constructions by default.
  if depth is None:
    depth = 5  # increased from 4 to probe stronger behavior

  # Pattern banks to introduce variety across recursion levels.
  START_SETS = [
    [2, 6, 10, 14],
    [3, 7, 11, 15],
    [4, 8, 12, 16],
    [5, 9, 13, 17],
  ]
  BRIDGE_SETS = [
    [(1, 5), (12, 16), (4, 9), (8, 13)],
    [(2, 6), (13, 17), (5, 10), (9, 14)],
    [(3, 7), (14, 18), (6, 11), (10, 15)],
    [(4, 8), (15, 19), (7, 12), (11, 16)],
  ]

  # Cyclic small perturbations of the scaling factor to break perfect self-similarity.
  SCALING_LIST = [1.00, 1.12, 0.95, 1.18, 1.03]

  T = [(0.0, 1.0)]
  for k in range(depth):
    lo = min(l for l, r in T)
    hi = max(r for l, r in T)
    base_delta = hi - lo if hi - lo > 0 else 1.0
    gamma = SCALING_LIST[k % len(SCALING_LIST)]
    delta = base_delta * gamma
    S = []
    # Use a rotating start set to diversify the expanded copies of T
    start_list = START_SETS[k % len(START_SETS)]
    for start in start_list:
      S += [(delta * start + l - lo, delta * start + r - lo) for l, r in T]
    # Use a rotating bridge pattern to insert cross-links without inflating omega
    bridge_set = BRIDGE_SETS[k % len(BRIDGE_SETS)]
    S += [(delta * a, delta * b) for (a, b) in bridge_set]
    T = S

    # Lightweight deterministic shrinker after certain levels:
    # if T has grown large in count, compress coordinates slightly to keep the span
    # from exploding while preserving relative overlaps.
    if len(T) > 800:
      # compress factor in (0.7, 0.9] depending on how large it is
      compress = 0.9 if len(T) < 1200 else 0.8
      # compute current lo/hi and compress towards center
      c_lo = min(l for l, r in T)
      c_hi = max(r for l, r in T)
      center = (c_lo + c_hi) / 2.0
      newT = []
      for (l, r) in T:
        nl = center + (l - center) * compress
        nr = center + (r - center) * compress
        # ensure non-degenerate
        if nr <= nl:
          nr = nl + 1.0
        newT.append((nl, nr))
      T = newT

  # Final mild normalization: shift to positive and round endpoints to integers
  # but preserve open-interval semantics by ensuring r > l.
  lo = min(l for l, r in T)
  if lo <= 0:
    shift = -lo + 1.0
  else:
    shift = 0.0
  normalized = []
  for (l, r) in T:
    L = int(round(l + shift))
    R = int(round(r + shift))
    if R <= L:
      R = L + 1
    normalized.append((L, R))
  return normalized

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