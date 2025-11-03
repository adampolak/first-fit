# EVOLVE-BLOCK-START

from math import gcd

def _normalize_grid(intervals):
  """
  Normalize endpoints to a compact integer grid while preserving order.
  Each unique endpoint is mapped to an increasing even integer.
  Returns a new list of (l, r) with integer coordinates.
  """
  endpoints = sorted(set([x for seg in intervals for x in seg]))
  coord = {}
  cur = 0
  for e in endpoints:
    coord[e] = cur
    cur += 2  # keep even gaps

  return [(coord[l], coord[r]) for (l, r) in intervals]


def construct_intervals(iterations=4, normalize=True):
  """
  Build a sequence of open intervals that forces FirstFit
  to use many colors while keeping the clique size controlled.
  This is a cross-over of:
    - Figure-4 style recursive expansion (four copies per level, plus four blockers)
    - a normalization step to map endpoints to a compact integer grid

  Arguments:
    iterations: number of recursive expansion steps (default 4)
    normalize: whether to apply endpoint normalization (default True)

  Returns:
    intervals: list of tuples (l, r) with integer endpoints (open intervals)
  """
  # Base gadget: a canonical unit interval
  T = [(0.0, 1.0)]

  # Recursive expansion: at each level replace T by copies plus blockers,
  # with an extra copy on the first iteration to heighten initial conflict.
  for i in range(iterations):
    lo = min(l for l, r in T)
    hi = max(r for l, r in T)
    delta = hi - lo
    S = []
    # Place four scaled/translated copies at staggered offsets
    for start in (2, 6, 10, 14):
      offset = delta * start - lo
      S.extend([(offset + l, offset + r) for l, r in T])
    # Add the four connecting long intervals (these create the coupling
    # between copies but are designed not to increase the clique beyond
    # the intended constant)
    S += [
      (delta * 1,  delta * 5),
      (delta * 12, delta * 16),
      (delta * 4,  delta * 9),
      (delta * 8,  delta * 13)
    ]
    T = S

  if normalize:
    return _normalize_grid(T)
  else:
    return T


# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()