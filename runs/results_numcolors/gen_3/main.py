# EVOLVE-BLOCK-START

def construct_intervals(depth=4, branching=5):
  """
  Parameterized recursive construction inspired by Kierstead-Trotter / Figure 4.

  Parameters:
    depth (int): recursion depth (default 4). Keep moderate to limit intervals.
    branching (int): number of scaled copies per level (default 5).

  Returns:
    list of (l, r) integer tuples representing open intervals in presentation order.
  """
  if depth <= 0:
    return [(0, 1)]

  T = [(0, 1)]
  for _ in range(depth):
    lo = min(l for l, r in T)
    hi = max(r for l, r in T)
    span = hi - lo
    S = []
    # starts: 2,6,10,... (spaced by 4) for branching copies
    starts = [2 + 4 * i for i in range(branching)]
    for s in starts:
      base = s * span - lo
      for (l, r) in T:
        S.append((int(l + base), int(r + base)))
    # adjacent connectors between blocks to propagate FF colors
    for i in range(len(starts) - 1):
      a = (1 + 4 * i) * span
      b = (5 + 4 * i) * span
      S.append((int(a), int(b)))
    # a couple of cross connectors (keeps clique small but couples blocks)
    if branching >= 3:
      a = (4 * 0 + 4) * span
      b = (4 * 2 + 9) * span
      S.append((int(a), int(b)))
    if branching >= 4:
      a = (4 * 1 + 8) * span
      b = (4 * 3 + 13) * span
      S.append((int(a), int(b)))
    T = S
  return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()