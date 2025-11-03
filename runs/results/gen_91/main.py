# EVOLVE-BLOCK-START

def construct_intervals(iterations: int = 4):
  """
  Construct a sequence of intervals of real line,
  in the order in which they are presented to FirstFit,
  so that it maximizes the number of colors used by FirstFit
  divided by the maximum number of intervals that cover a single point

  The initial implementation uses the construction from
  Figure 4 in https://arxiv.org/abs/1506.00192

  Returns:
    intervals: list of tuples, each tuple (l, r) represents an open interval from l to r
  """

  # --- Lightweight evaluators (deterministic) ---

  def _firstfit_colors(intervals):
    """Simulate FirstFit online coloring; open intervals allow reuse if l >= last_end."""
    last_end = []
    for (l, r) in intervals:
      placed = False
      for i, le in enumerate(last_end):
        if l >= le:
          last_end[i] = r
          placed = True
          break
      if not placed:
        last_end.append(r)
    return len(last_end)

  def _clique_number(intervals):
    """Sweep-line for open intervals; process right endpoints before left on ties."""
    events = []
    for (l, r) in intervals:
      if l < r:
        events.append((l, +1))
        events.append((r, -1))
    # right endpoints first at ties
    events.sort(key=lambda e: (e[0], 0 if e[1] == -1 else 1))
    cur = best = 0
    for _, t in events:
      cur += t
      if cur > best:
        best = cur
    return best

  # --- Builders for one recursive level ---

  def _make_copies(T, offsets, mode, lo, delta):
    """Create translated copies of T at the given offsets; combine by mode."""
    per_copy = []
    for start in offsets:
      offset = delta * start - lo
      per_copy.append([(l + offset, r + offset) for (l, r) in T])

    if mode == 'bycopy':
      S = []
      for lst in per_copy:
        S.extend(lst)
      return S
    else:  # 'roundrobin'
      m = len(T)
      S = []
      for i in range(m):
        for lst in per_copy:
          S.append(lst[i])
      return S

  def _blockers(delta, template='A'):
    """Connector intervals coupling the four copies."""
    if template == 'A':
      return [
        (delta * 1,  delta * 5),
        (delta * 12, delta * 16),
        (delta * 4,  delta * 9),
        (delta * 8,  delta * 13),
      ]
    else:  # 'B' (slightly shifted)
      return [
        (delta * 1.5, delta * 5.5),
        (delta * 11.5, delta * 15.5),
        (delta * 3.5,  delta * 8.5),
        (delta * 7.5,  delta * 12.5),
      ]

  def _build(k, offsets, order_mode, blk_pos, connector='A'):
    """k-level four-copy expansion with specified offsets/order/blocker placement."""
    T = [(0.0, 1.0)]
    for _ in range(k):
      lo = min(l for l, r in T)
      hi = max(r for l, r in T)
      delta = hi - lo
      copies = _make_copies(T, offsets, order_mode, lo, delta)
      blks = _blockers(delta, connector)
      if blk_pos == 'before':
        T = blks + copies
      elif blk_pos == 'after':
        T = copies + blks
      else:  # 'middle'
        mid = len(copies) // 2
        T = copies[:mid] + blks + copies[mid:]
    return T

  # --- Deterministic parameter sweep (small) ---

  offset_sets = [
    (2, 6, 10, 14),
    (1, 5, 9, 13),
    (3, 7, 11, 15),
    (0, 4, 8, 12),
  ]
  order_modes = ['bycopy', 'roundrobin']
  blocker_positions = ['after', 'before', 'middle']
  connector_templates = ['A', 'B']
  depths = [max(2, iterations), max(3, iterations + 1)]  # typically try 4 and 5

  best_T = None
  best_ratio = -1.0
  best_n = None
  best_cols = 0
  best_om = 0

  for k in depths:
    for offs in offset_sets:
      for omode in order_modes:
        for bpos in blocker_positions:
          for ctpl in connector_templates:
            T = _build(k, offs, omode, bpos, ctpl)
            om = _clique_number(T)
            if om <= 0:
              continue
            cols = _firstfit_colors(T)
            ratio = cols / om
            better = False
            if ratio > best_ratio + 1e-12:
              better = True
            elif abs(ratio - best_ratio) <= 1e-12:
              n = len(T)
              if best_n is None or n < best_n:
                better = True
              elif n == best_n and cols > best_cols:
                better = True
            if better:
              best_ratio = ratio
              best_T = T
              best_n = len(T)
              best_cols = cols
              best_om = om

  # Fallback to baseline pattern if sweep fails (should not happen)
  if best_T is None:
    T = [(0.0, 1.0)]
    for _ in range(iterations):
      lo = min(l for l, r in T)
      hi = max(r for l, r in T)
      delta = hi - lo
      S = []
      for start in (2, 6, 10, 14):
        S += [(delta * start + l - lo, delta * start + r - lo) for l, r in T]
      S += [
        (delta * 1,  delta * 5),
        (delta * 12, delta * 16),
        (delta * 4,  delta * 9),
        (delta * 8,  delta * 13),
      ]
      T = S
    best_T = T

  # Order-preserving normalization: remap endpoints to compact integers.
  # This preserves all overlap relations for open intervals.
  endpoints = sorted(set(x for seg in best_T for x in seg))
  mapping = {x: i for i, x in enumerate(endpoints)}
  T_norm = [(mapping[l], mapping[r]) for (l, r) in best_T]
  return T_norm

  # return [  # Figure 3, OPT=2, FF=4
  #   (2,3),
  #   (6,7),
  #   (10,11),
  #   (14,15),
  #   (1,5),
  #   (12,16),
  #   (4,9),
  #   (8,13),
  # ]

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()