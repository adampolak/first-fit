# EVOLVE-BLOCK-START

def construct_intervals(iterations=3):
  """
  Search over several deterministic variants of the 4-copy recursive gadget
  and return the best sequence of intervals (presented in arrival order)
  found within the search space. The returned intervals are normalized to
  integer endpoints (even integers) and are suitable for the existing
  evaluator that runs FirstFit on them.

  Arguments:
    iterations: nominal recursion depth (used as a hint; the search uses
                small depths {3,4} internally). This argument is kept for
                API compatibility.

  Returns:
    normalized: list of (l, r) integer tuples representing open intervals
                in the arrival order used to simulate FirstFit.
  """

  # Helper: build a recursive gadget given parameters.
  def build_gadget(base_gadget, k, starts_a, starts_b, conn_a, conn_b, mode='uniform'):
    T = [(float(l), float(r)) for (l, r) in base_gadget]
    for i in range(k):
      lo = min(l for l, r in T)
      hi = max(r for l, r in T)
      delta = hi - lo
      S = []
      # choose patterns depending on mode (uniform or alternate)
      if mode == 'uniform':
        starts = starts_a
        connectors = conn_a
      else:
        starts = starts_a if (i % 2 == 0) else starts_b
        connectors = conn_a if (i % 2 == 0) else conn_b
      # place copies
      for start in starts:
        offset = delta * start - lo
        for (l, r) in T:
          S.append((offset + l, offset + r))
      # add connector intervals (scaled by delta)
      for (a, b) in connectors:
        S.append((delta * a, delta * b))
      T = S
    return T

  # FirstFit online coloring simulation (open intervals).
  def first_fit_color(intervals):
    # last_end[color_index] = right endpoint of most recent interval assigned to that color
    last_end = []
    for l, r in intervals:
      placed = False
      for i in range(len(last_end)):
        if last_end[i] <= l:  # open intervals: equal endpoint is allowed
          last_end[i] = r
          placed = True
          break
      if not placed:
        last_end.append(r)
    return len(last_end)

  # Compute clique number (maximum number intervals covering any point).
  def clique_number(intervals):
    # For open intervals, an endpoint is not included. Process end events before start events.
    events = []
    for l, r in intervals:
      events.append((l, +1))   # start
      events.append((r, -1))   # end
    # sort by coordinate; for tie: end (-1) should come before start (+1) => natural ordering suffices
    events.sort(key=lambda x: (x[0], x[1]))
    cur = 0
    best = 0
    for x, t in events:
      if t == -1:
        cur -= 1
      else:
        cur += 1
        if cur > best:
          best = cur
    return best

  # Normalize to increasing even integers preserving order and distinctness.
  def normalize_intervals(intervals):
    eps = 1e-12
    endpoints = sorted(set([x for seg in intervals for x in seg]))
    coord = {}
    cur = 0
    for e in endpoints:
      coord[e] = cur
      cur += 2
    normalized = [(coord[l], coord[r]) for (l, r) in intervals]
    return normalized

  # Candidate parameter libraries (deterministic)
  base_gadgets = [
    [(0.0, 1.0)],                  # canonical base
    [(0.0, 1.0), (3.0, 4.0)]       # slightly richer base gadget
  ]

  starts_options = [
    (2, 6, 10, 14),
    (1, 5, 9, 13),
    (3, 7, 11, 15),
    (0, 4, 8, 12)
  ]

  connector_templates = [
    # original Figure-4 connectors
    [(1.0, 5.0), (12.0, 16.0), (4.0, 9.0), (8.0, 13.0)],
    # shifted by 0.5 to break alignment
    [(1.5, 5.5), (12.5, 16.5), (4.5, 9.5), (8.5, 13.5)],
    # alternative templates (from recommendations)
    [(0.0, 4.0), (11.0, 15.0), (3.0, 8.0), (7.0, 12.0)],
    [(0.0, 6.0), (10.0, 16.0), (5.0, 11.0), (9.0, 15.0)]
  ]

  # Search parameters (kept modest so construction is small and evaluation fast)
  k_values = [3, 4]     # try depth 3 and 4
  modes = ['uniform', 'alternate']

  best_candidate = None
  best_metric = -1.0
  best_info = None

  # Enumerate deterministic candidate space
  for base in base_gadgets:
    for k in k_values:
      for si in range(len(starts_options)):
        starts_a = starts_options[si]
        starts_b = starts_options[(si + 1) % len(starts_options)]
        for ci in range(len(connector_templates)):
          conn_a = connector_templates[ci]
          conn_b = connector_templates[(ci + 1) % len(connector_templates)]
          for mode in modes:
            # Build candidate gadget
            cand = build_gadget(base, k, starts_a, starts_b, conn_a, conn_b, mode=mode)
            # Normalize coordinates (but keep floats for coloring/clique computation)
            # We won't normalize yet because clique/FF depend only on order, but to keep
            # consistent evaluation we normalize at the end for final output.
            # Here compute metrics on raw floats.
            alg_cols = first_fit_color(cand)
            opt = clique_number(cand)
            # avoid divide by zero
            if opt <= 0:
              ratio = 0.0
            else:
              ratio = alg_cols / float(opt)
            # choose by best ratio, tie-breaker fewer intervals then smaller alg_cols
            if ratio > best_metric or (abs(ratio - best_metric) < 1e-12 and (best_candidate is None or (len(cand) < len(best_candidate) or (len(cand) == len(best_candidate) and alg_cols > (best_info or (0,0))[0])))):
              best_candidate = cand
              best_metric = ratio
              best_info = (alg_cols, opt, len(cand), (base, k, starts_a, conn_a, mode))

  # If search found nothing unusual, fall back to simple default construction
  if best_candidate is None:
    best_candidate = build_gadget([(0.0, 1.0)], iterations, starts_options[0], starts_options[1], connector_templates[0], connector_templates[1], mode='uniform')

  # Conservative pruning pass: try removing intervals one-by-one (reverse order) while
  # preserving or improving the FirstFit/OPT ratio
  def prune_intervals(intervals):
    cur_intervals = list(intervals)
    cur_alg = first_fit_color(cur_intervals)
    cur_opt = clique_number(cur_intervals)
    cur_ratio = cur_alg / max(1, cur_opt)
    # iterate in reverse to try removing later (often redundant) intervals first
    i = len(cur_intervals) - 1
    while i >= 0:
      candidate = cur_intervals[:i] + cur_intervals[i+1:]
      new_alg = first_fit_color(candidate)
      new_opt = clique_number(candidate)
      new_ratio = new_alg / max(1, new_opt)
      # Accept removal if ratio doesn't decrease (allow small floating error)
      if new_ratio + 1e-12 >= cur_ratio:
        cur_intervals = candidate
        cur_alg, cur_opt, cur_ratio = new_alg, new_opt, new_ratio
        # keep same index (next element shifts into current position)
      else:
        i -= 1
    return cur_intervals, cur_alg, cur_opt

  pruned, final_alg, final_opt = prune_intervals(best_candidate)

  # Normalize endpoints to even integers for tidy output
  normalized = normalize_intervals(pruned)

  return normalized

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()