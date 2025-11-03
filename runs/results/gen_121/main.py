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


def _compute_four_copies_and_blockers(T):
  """
  Given a current template T (list of intervals), produce the next level
  by placing four translated copies and appending four long blockers.
  Offsets used are canonical: 2, 6, 10, 14 (relative to current delta/lo).
  """
  # current span
  lo = min(l for l, r in T)
  hi = max(r for l, r in T)
  delta = hi - lo

  # four copies at offsets 2,6,10,14
  S = []
  for start in (2, 6, 10, 14):
    offset = delta * start - lo
    S.extend([(offset + l, offset + r) for (l, r) in T])

  # four connecting long intervals
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

  # Recursive expansion: at each level replace T by four copies plus blockers
  for _ in range(iterations):
    T = _compute_four_copies_and_blockers(T)

  # Inject a wave of short intervals to consume extra FirstFit colors without raising omega
  lo = min(l for l, r in T)
  hi = max(r for l, r in T)
  delta = hi - lo
  small_len = delta * 0.1
  for start in (2, 6, 10, 14):
    offset = delta * start - lo
    # place a short interval inside each copy region
    T.append((offset + 0.1 * delta, offset + 0.1 * delta + small_len))

  if normalize:
    return _normalize_grid(T)
  else:
    return T


# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()