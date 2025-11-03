# EVOLVE-BLOCK-START

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


def _compute_four_copies_and_blockers(T, starts):
  """
  Given a current template T (list of intervals), produce the next level
  by placing translated copies according to 'starts' offsets and appending four long blockers.
  starts: tuple of offsets (multipliers) for copies
  """
  lo = min(l for l, r in T)
  hi = max(r for l, r in T)
  delta = hi - lo

  S = []
  for start in starts:
    offset = delta * start - lo
    S.extend([(offset + l, offset + r) for (l, r) in T])

  S += [
    (delta * 1,  delta * 5),
    (delta * 12, delta * 16),
    (delta * 4,  delta * 9),
    (delta * 8,  delta * 13)
  ]
  return S


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

  # Recursive expansion: at each level replace T by four copies plus blockers with cycling patterns
  patterns = [(2, 6, 10, 14), (1, 5, 9, 13), (3, 7, 11, 15), (0, 4, 8, 12)]
  for i in range(iterations):
    offsets = patterns[i % len(patterns)]
    T = _compute_four_copies_and_blockers(T, offsets)

  if normalize:
    return _normalize_grid(T)
  else:
    return T


# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()