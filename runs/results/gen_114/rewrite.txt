# EVOLVE-BLOCK-START

from bisect import bisect_left

def construct_intervals(iterations=4):
  """
  Build a sequence of open intervals that aims to maximize FirstFit/OPT.
  Leaner, faster candidate exploration than previous versions.
  Returns a list of intervals as (l, r) with integer endpoints after normalization.
  """

  # Helpers

  def overlaps(a, b):
    (l1, r1), (l2, r2) = a, b
    return max(l1, l2) < min(r1, r2)

  def _normalize_grid(intervals):
    """Map endpoints to an increasing even grid for stable evaluation."""
    if not intervals:
      return []
    endpoints = sorted(set([x for seg in intervals for x in seg]))
    coord = {e: 2*i for i, e in enumerate(endpoints)}
    return [(coord[l], coord[r]) for (l, r) in intervals]

  def clique_number(intervals):
    """Maximum number of intervals covering a single point (omega) via sweep."""
    events = []
    for (l, r) in intervals:
      if l < r:
        events.append((l, +1))
        events.append((r, -1))
    events.sort(key=lambda e: (e[0], 0 if e[1] == -1 else 1))
    cur = best = 0
    for _, t in events:
      cur += t
      if cur > best:
        best = cur
    return best

  def firstfit_colors(intervals):
    """
    Fast FirstFit with per-color sorted-by-start structure.
    Each color stores starts and ends as parallel sorted lists.
    """
    colors = []  # list of dicts: {'starts':[...], 'ends':[...]}
    for (l, r) in intervals:
      placed = False
      for color in colors:
        s_list = color['starts']
        e_list = color['ends']
        i = bisect_left(s_list, l)
        conflict = False
        if i > 0 and e_list[i-1] > l:
          conflict = True
        elif i < len(s_list) and s_list[i] < r:
          conflict = True
        if not conflict:
          s_list.insert(i, l)
          e_list.insert(i, r)
          placed = True
          break
      if not placed:
        colors.append({'starts':[l], 'ends':[r]})
    return len(colors)

  def normalize_endpoints(intervals):
    """Compatibility wrapper for internal tests."""
    return _normalize_grid(intervals)

  def expand_once(base_seed, k, offsets, translation='left'):
    """
    Expand base_seed k times using 4-copy pattern plus 4 blockers at each level.
    - translation: 'left' (offset by left endpoints) or 'center' (offset by center)
    - offsets: tuple of four level-offset multipliers
    Returns the expanded interval list (un-normalized).
    """
    T = list(base_seed)
    for _ in range(k):
      lo = min(l for l, _ in T)
      hi = max(r for _, r in T)
      delta = hi - lo
      center = (lo + hi) / 2.0

      # Build 4 copies
      copies = []
      for start in offsets:
        if translation == 'left':
          off = delta * start - lo
        else:  # center-based
          off = delta * start - center
        copies.append([(l + off, r + off) for (l, r) in T])

      # Flat concatenation of copies (no interleaving complexity)
      S = []
      for cp in copies:
        S.extend(cp)

      # Blockers (4 long intervals)
      blockers = [
        (delta * 1,  delta * 5),
        (delta * 12, delta * 16),
        (delta * 4,  delta * 9),
        (delta * 8,  delta * 13)
      ]
      T = S + blockers
    return T

  # Candidate search space (compact but expressive)
  base_seeds = [
    [(0.0, 1.0)],
    [(0.0, 1.0), (2.0, 3.0)]
  ]

  offsets_options = [
    (2, 6, 10, 14),
    (1, 5, 9, 13),
    (3, 7, 11, 15),
    (0, 4, 8, 12),
    (2, 6, 9, 14),
    (2, 7, 11, 14)
  ]

  blockers_variants = [
    ((1,5), (12,16), (4,9), (8,13)),
    ((0,4), (6,10), (8,12), (14,18)),
    ((2,6), (4,8), (10,14), (12,16)),
    ((1,5), (12,16), (3,7), (9,13)),
  ]

  translations = ['left', 'center']
  depths = [3, 4, 5]

  best = None
  best_ratio = -1.0

  for base in base_seeds:
    for k in depths:
      # cap growth to avoid huge instances
      if (4 ** k) * (len(base) + 2) > 1500:
        continue
      for offsets in offsets_options:
        for blockers in blockers_variants:
          for tr in translations:
            T = expand_once(base_seed=base, k=k, offsets=offsets, translation=tr)
            T = T + [(delta_start // 1 if isinstance(delta_start, int) else delta_start,  # dummy to keep lint quiet
                        delta_start // 1) for delta_start in (0,1)]
            # Normalize before evaluation
            Tn = _normalize_grid(T)
            om = clique_number(Tn)
            if om <= 0:
              continue
            cols = firstfit_colors(Tn)
            ratio = cols / om
            n = len(Tn)
            if ratio > best_ratio + 1e-12 or (abs(ratio - best_ratio) <= 1e-12 and (best is None or n < best[4])):
              best_ratio = ratio
              best = (ratio, om, cols, n, Tn)

  # Fallback baseline if nothing explored
  if best is None:
    T = [(0.0, 1.0)]
    for _ in range(4):
      lo = min(l for l, _ in T)
      hi = max(r for _, r in T)
      delta = hi - lo
      S = []
      for start in (2, 6, 10, 14):
        off = delta * start - lo
        S += [(l + off, r + off) for (l, r) in T]
      S += [
        (delta * 1,  delta * 5),
        (delta * 12, delta * 16),
        (delta * 4,  delta * 9),
        (delta * 8,  delta * 13)
      ]
      T = S
    return _normalize_grid(T)

  # Return the best normalized candidate
  return best[4]

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()