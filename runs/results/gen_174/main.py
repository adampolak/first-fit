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


def _compute_four_copies_and_blockers(T, extra_offset=False, alt_offsets=False):
  """
  Given a current template T (list of intervals), produce the next level
  by placing translated copies and appending long blockers.
  - alt_offsets=True switches base offsets to (1,5,9,13)
  - extra_offset=True adds a fifth copy at multiplier 18
  """
  # current span
  lo = min(l for l, r in T)
  hi = max(r for l, r in T)
  delta = hi - lo

  # determine base offsets
  if alt_offsets:
    base = (1, 5, 9, 13)
  else:
    base = (2, 6, 10, 14)
  offsets = base + ((18,) if extra_offset else ())

  # copies at specified offsets
  S = []
  for start in offsets:
    offset = delta * start - lo
    for l, r in T:
      S.append((offset + l, offset + r))

  # four connecting long intervals (blockers)
  S += [
    (delta * 1,  delta * 5),
    (delta * 12, delta * 16),
    (delta * 4,  delta * 9),
    (delta * 8,  delta * 13),
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

  # Recursive expansion: alternate offset patterns and optionally add an extra copy on the first iteration
  for i in range(iterations):
    extra = (i == 0)
    alt = (i % 2 == 1)
    T = _compute_four_copies_and_blockers(T, extra_offset=extra, alt_offsets=alt)

  if normalize:
    return _normalize_grid(T)
  else:
    return T


# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()