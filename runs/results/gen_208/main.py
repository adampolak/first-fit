# EVOLVE-BLOCK-START

def construct_intervals(iterations=4):
  """
  Construct a sequence of open intervals presented to FirstFit that reproduces
  the canonical Figure-4 adversary (four copies plus four long connectors per level),
  while eliminating the extra first-iteration copy and normalizing endpoints onto
  a compact integer grid. This preserves omega=5 and the known 13-color FirstFit
  behavior, but reduces the instance size and removes floating artifacts.
  """

  # Base gadget
  T = [(0.0, 1.0)]

  # Each iteration: four translated/scaled copies plus four connectors
  for _ in range(iterations):
    lo = min(l for l, r in T)
    hi = max(r for l, r in T)
    delta = hi - lo
    S = []
    # Four tiled copies at offsets 2, 6, 10, 14 (canonical placement)
    for start in (2, 6, 10, 14):
      offset = delta * start - lo
      S.extend([(offset + l, offset + r) for (l, r) in T])
    # Four long connectors (blockers) as in Figure 4
    S += [
      (delta * 1,  delta * 5),
      (delta * 12, delta * 16),
      (delta * 4,  delta * 9),
      (delta * 8,  delta * 13),
    ]
    T = S

  # Normalize endpoints to a compact even-integer grid
  endpoints = sorted(set([x for seg in T for x in seg]))
  coord = {}
  cur = 0
  for e in endpoints:
    coord[e] = cur
    cur += 2
  return [(coord[l], coord[r]) for (l, r) in T]

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()