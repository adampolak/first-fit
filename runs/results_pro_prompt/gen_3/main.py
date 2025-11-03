# EVOLVE-BLOCK-START

def construct_intervals():
  """
  Construct a sequence of open intervals presented to FirstFit.
  Improved version: search over recursion depths and arrival orders
  of a fractal/tile construction to maximize FirstFit / OPT ratio.

  Returns:
    intervals: list of tuples (l, r) representing open intervals [l, r)
  """
  import random

  def build_fractal(depth, starts=(2, 6, 10, 14), caps=((1, 5), (12, 16), (4, 9), (8, 13))):
    T = [(0.0, 1.0)]
    for _ in range(depth):
      lo = min(l for l, r in T)
      hi = max(r for l, r in T)
      delta = hi - lo
      S = []
      for start in starts:
        shift = delta * start - lo
        for l, r in T:
          S.append((shift + l, shift + r))
      for a, b in caps:
        S.append((delta * a, delta * b))
      T = S
    # Normalize so leftmost point is 0.0 (keeps numbers small)
    lo = min(l for l, r in T)
    T = [(l - lo, r - lo) for l, r in T]
    return T

  def firstfit_for_order(intervals):
    # Greedy FirstFit in the given order
    colors = []
    for i, (l, r) in enumerate(intervals):
      used = set()
      # Check overlap with previously colored intervals
      # For open intervals (l, r) overlaps (l2, r2) iff l < r2 and l2 < r
      for j, (lj, rj) in enumerate(intervals[:i]):
        if lj < r and l < rj:
          used.add(colors[j])
      # assign smallest positive integer not used
      c = 1
      while c in used:
        c += 1
      colors.append(c)
    return colors

  def firstfit_colors(intervals):
    assigned = firstfit_for_order(intervals)
    return max(assigned) if assigned else 0

  def clique_number(intervals):
    # Sweep-line to compute maximum number of intervals covering a point
    # For open intervals, endpoints do not count as overlap; process ends before starts on ties.
    events = []
    for l, r in intervals:
      events.append((l, 1))   # start
      events.append((r, -1))  # end
    events.sort(key=lambda x: (x[0], x[1]))  # end (-1) sorted before start (1) on ties
    curr = 0
    best = 0
    for pos, delta in events:
      curr += delta
      if curr > best:
        best = curr
    return best

  # Try several depths and several arrival orderings; pick best ratio
  best_seq = None
  best_meta = None
  best_ratio = -1.0

  random_seed = 1234567

  # Expand depth search to explore more patterns (depth up to 5)
  for depth in range(1, 6):
    seqs = [
      ("fractal_A", build_fractal(depth)),
      ("fractal_B", build_fractal(depth, starts=(1, 5, 9, 13)))
    ]
    for base_name, seq in seqs:
      opt = clique_number(seq) or 1

      # Generate a diverse collection of arrival orders (more seeds and center-based orders)
      orders = []
      orders.append((f"{base_name}", list(seq)))
      orders.append(("reversed", list(reversed(seq))))
      orders.append(("short_first", sorted(seq, key=lambda x: (x[1] - x[0], x[0]))))
      orders.append(("long_first", sorted(seq, key=lambda x: -(x[1] - x[0]))))
      orders.append(("left_first", sorted(seq, key=lambda x: x[0])))
      orders.append(("right_first", sorted(seq, key=lambda x: -x[0])))
      orders.append(("center", sorted(seq, key=lambda x: (x[0] + x[1]) / 2.0)))
      orders.append(("center_rev", sorted(seq, key=lambda x: -(x[0] + x[1]) / 2.0)))
      # Deterministic random seeds for diversity
      rng = random.Random(random_seed + depth * 10)
      for i in range(6):
        tmp = list(seq)
        rng.shuffle(tmp)
        orders.append((f"rand_{base_name}_{i}", tmp))

      for name, candidate in orders:
        colors_used = firstfit_colors(candidate)
        ratio = colors_used / opt
        if ratio > best_ratio + 1e-12:
          best_ratio = ratio
          best_seq = candidate
          best_meta = (depth, name, colors_used, opt, len(candidate))

  # Fallback (shouldn't happen): return depth=2 original pattern
  if best_seq is None:
    best_seq = build_fractal(2)
    best_meta = (2, "fallback", firstfit_colors(best_seq), clique_number(best_seq), len(best_seq))

  # Print a short validation line
  depth, order_name, ff_colors, opt_val, n_intervals = best_meta
  print(f"construct_intervals: depth={depth}, order={order_name}, n={n_intervals}, FirstFit={ff_colors}, OPT={opt_val}, ratio={ff_colors/opt_val:.3f}")

  return best_seq

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()