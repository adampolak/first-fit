# EVOLVE-BLOCK-START

def construct_intervals(iterations=4):
  """
  Construct a sequence of intervals of real line to maximize FirstFit/omega.
  We explore small variations inspired by adversarial arrival-order engineering:
  - schedule of blockers vs copies at each level (before/after/split)
  - optional extra copy at the first iteration to couple layers
  We evaluate candidates via an in-block FirstFit and omega computation, and return the best.

  Enhancements in this edit:
   * shrink_prune: conservative one-by-one removal of intervals (prefer long ones)
     while preserving the observed FirstFit/OPT ratio. This reduces interval count n
     and improves the combined score without weakening the witness.
   * _normalize_grid: map unique endpoints to a compact integer grid (even integers).
  """

  # --- Helper routines (kept lightweight and deterministic) ---

  def firstfit_colors(intervals):
    """Simulate FirstFit using per-color last end-point tracking."""
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

  def clique_number(intervals):
    """Sweep for open intervals; process right endpoints before left at ties."""
    events = []
    for (l, r) in intervals:
      if l < r:
        events.append((l, +1))
        events.append((r, -1))
    # ensure right-endpoints decrease before left at same coordinate
    events.sort(key=lambda e: (e[0], 0 if e[1] == -1 else 1))
    cur = best = 0
    for _, t in events:
      cur += t
      if cur > best:
        best = cur
    return best

  def _normalize_grid(intervals):
    """Map unique endpoints to increasing even integers preserving order."""
    endpoints = sorted(set([x for seg in intervals for x in seg]))
    coord = {}
    cur = 0
    for e in endpoints:
      coord[e] = cur
      cur += 2
    return [(coord[l], coord[r]) for (l, r) in intervals]

  def shrink_prune(intervals, target_ratio):
    """
    Try to remove intervals one-by-one (prefer longer intervals first).
    Accept a removal only if the FirstFit/OPT ratio remains >= target_ratio.
    Conservative: stops when no single removal is safe.
    """
    cur = list(intervals)
    # order candidates by decreasing length (prefer removing long blockers)
    def lengths_order(lst):
      return [i for i, _ in sorted(((idx, seg[1]-seg[0]) for idx, seg in enumerate(lst)),
                                   key=lambda x: (-x[1], x[0]))]
    order = lengths_order(cur)
    changed = True
    while changed:
      changed = False
      for idx in order:
        if idx >= len(cur):
          continue
        cand = cur[:idx] + cur[idx+1:]
        if not cand:
          continue
        alg = firstfit_colors(cand)
        opt = clique_number(cand)
        if opt == 0:
          continue
        ratio = alg / opt
        if ratio >= target_ratio - 1e-12:
          # accept removal
          cur = cand
          order = lengths_order(cur)
          changed = True
          break
    return cur

  def build_candidate(k, schedule='after', extra_first=False):
    """
    Build k-level recursive pattern using four scaled copies per level,
    with four long connectors ("blockers"), and schedule controlling
    the arrival order within each level.
    schedule in {'after','before','split'}.
    """
    T = [(0.0, 1.0)]
    for i in range(k):
      lo = min(l for l, r in T)
      hi = max(r for l, r in T)
      delta = hi - lo

      # Offsets for this level; optionally add a fifth copy on the first iteration
      offs = (2, 6, 10, 14, 18) if (extra_first and i == 0) else (2, 6, 10, 14)

      # Prepare copies
      def make_copies(from_T, offsets):
        S = []
        for start in offsets:
          offset = delta * start - lo  # left-anchored translation
          for (l, r) in from_T:
            S.append((l + offset, r + offset))
        return S

      # Blockers as in Figure 4
      blockers = [
        (delta * 1,  delta * 5),
        (delta * 12, delta * 16),
        (delta * 4,  delta * 9),
        (delta * 8,  delta * 13),
      ]

      if schedule == 'after':
        S = make_copies(T, offs) + blockers
      elif schedule == 'before':
        S = list(blockers) + make_copies(T, offs)
      else:  # 'split': half the copies, then blockers, then remaining copies
        h = len(offs) // 2
        first = offs[:h]
        second = offs[h:]
        S = make_copies(T, first) + list(blockers) + make_copies(T, second)

      T = S
    return T

  # --- Enumerate a modest set of candidates and select the best by ratio ---
  depths = [max(2, iterations - 1), iterations]  # try one less and the given
  schedules = ['after', 'before', 'split']
  extras = [False, True]

  best_T = None
  best_ratio = -1.0
  best_n = None
  best_cols = 0
  best_om = 0

  for k in depths:
    for sch in schedules:
      for ex in extras:
        T = build_candidate(k, schedule=sch, extra_first=ex)
        om = clique_number(T)
        if om == 0:
          continue
        cols = firstfit_colors(T)
        ratio = cols / om
        n = len(T)
        # Prefer higher ratio; tie-break by fewer intervals then more colors
        better = False
        if ratio > best_ratio + 1e-12:
          better = True
        elif abs(ratio - best_ratio) <= 1e-12:
          if best_n is None or n < best_n:
            better = True
          elif n == best_n and cols > best_cols:
            better = True
        if better:
          best_ratio = ratio
          best_T = T
          best_n = n
          best_cols = cols
          best_om = om

  # Fallback to baseline 4-iteration construction if search fails
  if best_T is None:
    T = [(0.0, 1.0)]
    for i in range(iterations):
      lo = min(l for l, r in T)
      hi = max(r for l, r in T)
      delta = hi - lo
      S = []
      for start in (2, 6, 10, 14):
        S += [(delta * start + l - lo, delta * start + r - lo) for l, r in T]
      S += [
        (delta * 1, delta * 5),
        (delta * 12, delta * 16),
        (delta * 4, delta * 9),
        (delta * 8, delta * 13)
      ]
      T = S
    best_T = T
    best_ratio = firstfit_colors(best_T) / max(1, clique_number(best_T))

  # Apply conservative pruning to reduce n while preserving ratio
  pruned = shrink_prune(best_T, best_ratio)

  # Normalize endpoints to compact integer grid for final output
  final = _normalize_grid(pruned)

  return final

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