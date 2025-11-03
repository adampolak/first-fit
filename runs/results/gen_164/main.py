# EVOLVE-BLOCK-START

def construct_intervals(iterations: int = 4):
  """
  Construct a sequence of intervals of real line,
  in the order in which they are presented to FirstFit,
  so that it maximizes the number of colors used by FirstFit
  divided by the maximum number of intervals that cover a single point

  The variant cycles through four canonical copy patterns across iterations
  to increase adversarial diversity while preserving omega.
  Patterns (per-iteration) are:
  A: (2, 6, 10, 14)
  B: (1, 5, 9, 13)
  C: (3, 7, 11, 15)
  D: (0, 4, 8, 12)
  The cycle repeats every four iterations.

  Returns:
    intervals: list of tuples, each tuple (l, r) represents an open interval from l to r
  """

  # Classic 4-copy/4-blocker recursive pattern (baseline adversary) with cycling offsets
  T = [(0.0, 1.0)]
  patterns = [
    (2, 6, 10, 14),  # A
    (1, 5, 9, 13),   # B
    (3, 7, 11, 15),  # C
    (0, 4, 8, 12),   # D
  ]
  for i in range(iterations):
    lo = min(l for l, r in T)
    hi = max(r for l, r in T)
    delta = hi - lo
    S = []
    pat = patterns[i % 4]
    for start in pat:
      S += [(delta * start + l - lo, delta * start + r - lo) for l, r in T]
    S += [
      (delta * 1,  delta * 5),
      (delta * 12, delta * 16),
      (delta * 4,  delta * 9),
      (delta * 8,  delta * 13),
    ]
    T = S

  # Order-preserving normalization: remap endpoints to compact integers.
  endpoints = sorted(set(x for seg in T for x in seg))
  mapping = {x: i for i, x in enumerate(endpoints)}
  T_norm = [(mapping[l], mapping[r]) for (l, r) in T]
  return T_norm

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()