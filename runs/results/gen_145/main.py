# EVOLVE-BLOCK-START

def construct_intervals(iterations=4):
  """
  Construct a sequence of open intervals (l, r) presented in order to FirstFit.
  The construction is a recursive expansion inspired by Figure 4 of
  https://arxiv.org/abs/1506.00192. Increasing 'iterations' strengthens the
  adversary and typically increases the FirstFit/OPT ratio.

  Arguments:
    iterations: number of recursive expansion steps (default 4)

  Returns:
    intervals: list of tuples (l, r) with integer endpoints (open intervals)
  """

  # Base gadget: a canonical small pattern (one unit interval)
  T = [(0.0, 1.0)]

  # For each iteration we replace T by four translated/scaled copies
  # plus four long connector intervals (as in the Figure-4 construction).
  # This preserves a small clique size while amplifying FirstFit colors.
  for _ in range(iterations):
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

  # Add wave intervals to boost FirstFit colors without increasing clique number
  lo = min(l for l, r in T)
  hi = max(r for l, r in T)
  span = hi - lo
  # choose wave length as small fraction of total span
  wave_length = span / 100.0
  # determine number of non‐overlapping waves, capped for compactness
  wave_count = min(int(span / (wave_length * 3.0)), 60)
  for i in range(wave_count):
    start = lo + (i + 1) * span / (wave_count + 1)
    T.append((start, start + wave_length))
  # Normalize endpoints to a compact integer grid while preserving order.
  # Map each unique endpoint to an increasing even integer (to ensure
  # positive lengths and avoid degeneracy).
  eps = 1e-9
  endpoints = sorted(set([x for seg in T for x in seg]))
  # Use a dictionary mapping each float endpoint to an integer coordinate.
  coord = {}
  cur = 0
  for e in endpoints:
    # assign increasing even integers (spacing by 2)
    coord[e] = cur
    cur += 2

  normalized = [(coord[l], coord[r]) for (l, r) in T]

  return normalized

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()