# EVOLVE-BLOCK-START

def construct_intervals(iterations=4):
  """
  Construct a sequence of intervals of real line to maximize FirstFit/omega.
  Improvements over prior implementations:
   - deterministic exploration over left/center translations and multiple offset patterns
   - normalization to a compact integer grid
   - deterministic pruning to shrink witness while preserving the score
  """

  # --- Helpers (deterministic, lightweight) ---

  def _normalize_grid(intervals):
    """
    Map all unique endpoints to increasing even integers (0,2,4,...).
    Keeps numbers small and integral, ensuring stable scoring.
    """
    if not intervals:
      return []
    endpoints = sorted(set(x for seg in intervals for x in seg))
    coord = {}
    cur = 0
    for e in endpoints:
      coord[e] = cur
      cur += 2
    return [(coord[l], coord[r]) for (l, r) in intervals]

  def overlaps(a, b):
    (l1, r1), (l2, r2) = a, b
    return max(l1, l2) < min(r1, r2)

  def firstfit_colors(intervals):
    """
    Efficient FirstFit using per-color last-end tracking.
    Works because in a valid FirstFit coloring, intervals assigned to a color
    are non-overlapping.
    """
    last_end = []  # last_end[c-1] = right endpoint of last interval in color c
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
    """
    Open-interval omega via sweep-line. Endpoints handled so that open-interval overlaps counted.
    """
    events = []
    for (l, r) in intervals:
      if l < r:
        events.append((l, +1))
        events.append((r, -1))
    events.sort(key=lambda e: (e[0], 0 if e[1] == -1 else 1))
    cur = 0
    best = 0
    for _, t in events:
      cur += t
      if cur > best:
        best = cur
    return best

  def build_pattern(k, base_seed, offsets, blockers, translation, blocker_anchor, extra_copies=0):
    """
    Recursively expand the base_seed k times using the 4-copy + 4-blocker scheme.
    translation: 'left' or 'center' controls how copies are positioned.
    offsets: tuple of integer multipliers for the level translation.
    blockers: 4-blocker geometry as 4 (a,b) pairs in each iteration.
    extra_copies: optionally add an extra copy on the first level to diversify layout.
    """
    T = list(base_seed)
    for i in range(k):
      lo = min(l for l, r in T)
      hi = max(r for l, r in T)
      delta = hi - lo
      center = (lo + hi) / 2.0

      offs = tuple(list(offsets) + ([18] if (extra_copies and i == 0) else []))

      # helper: place copies
      S = []
      for start in offs:
        if translation == 'left':
          offset = delta * start - lo
        else:  # center-based
          offset = delta * start - center
        for (l, r) in T:
          S.append((l + offset, r + offset))

      # blockers: four long intervals that tie the copies together
      blockers_this = [
        (delta * 1,  delta * 5),
        (delta * 12, delta * 16),
        (delta * 4,  delta * 9),
        (delta * 8,  delta * 13),
      ]

      S = S + blockers_this

      T = S

    return T

  def evaluate(intervals):
    """
    Compute normalized witness, omega, FirstFit colors and ratio (cols/omega).
    Returns (score, omega, colors, n, normalized_intervals)
    """
    Tn = _normalize_grid(intervals)
    n = len(Tn)
    if n == 0:
      return (-1.0, 0, 0, n, Tn)
    om = clique_number(Tn)
    if om == 0:
      return (-1.0, 0, 0, n, Tn)
    cols = firstfit_colors(Tn)
    ratio = cols / om
    # small penalty for larger witnesses to prefer compact witnesses when ratio ties
    score = ratio - 1e-6 * (n / 10000.0)
    return (score, om, cols, n, Tn)

  def prune_by_ratio(intervals, target_ratio):
    """
    Deterministic pruning:
    greedily remove intervals (long blockers first) while preserving ratio >= target_ratio.
    """
    cur = list(intervals)
    if not cur:
      return cur
    improved = True
    # order by decreasing length
    while improved:
      improved = False
      order = sorted(range(len(cur)), key=lambda i: (-(cur[i][1] - cur[i][0]), i))
      for idx in order:
        cand = cur[:idx] + cur[idx+1:]
        if not cand:
          continue
        s, om, cols, n, _ = evaluate(cand)
        if om > 0:
          ratio = cols / om
          if ratio >= target_ratio - 1e-12:
            cur = cand
            improved = True
            break
    return cur

  # --- Enumerate a modest but meaningful blueprint space deterministically ---

  offsets_set = [
    (2, 6, 10, 14),  # baseline
    (1, 5, 9, 13),   # shifted pattern
    (3, 7, 11, 15),  # alternate stagger
  ]
  blockers_templates = [
    ((1, 5), (12, 16), (4, 9), (8, 13)),  # Template A
    ((0, 4), (6, 10), (8, 12), (14, 18)), # Template B
    ((2, 6), (4, 8), (10, 14), (12, 16)), # Template C
  ]
  translations = ['left', 'center']
  blocker_anchors = ['left', 'center']
  depths = [3, 4]  # keep search light but informative
  base_seeds = [
    [(0.0, 1.0)],
    [(0.0, 1.0), (2.0, 3.0)]
  ]
  # small extra copy option to diversify layout
  extra_options = [0, 1]

  best = None  # (score, omega, colors, n, normalized_intervals, raw_intervals)

  for base in base_seeds:
    for k in depths:
      # skip overly large combos by construction size estimate
      if (4 ** k) * (len(base) + 2) > 2000:
        continue
      for offsets in offsets_set:
        for blockers in blockers_templates:
          for translation in translations:
            for anchor in blocker_anchors:
              for extra in extra_options:
                T = build_pattern(
                    k=k,
                    base_seed=base,
                    offsets=offsets,
                    blockers=blockers,
                    translation=translation,
                    blocker_anchor=anchor,
                    extra_copies=extra
                )
                score, om, cols, n, Tn = evaluate(T)
                cand = (score, om, cols, n, Tn, T)
                if best is None:
                  best = cand
                else:
                  if cand[0] > best[0] + 1e-12:
                    best = cand
                  elif abs(cand[0] - best[0]) <= 1e-12:
                    if cand[3] < best[3]:
                      best = cand
                    elif cand[3] == best[3] and cand[2] > best[2]:
                      best = cand

  # Fallback to a known baseline if search failed
  if best is None:
    k = 4
    T = [(0.0, 1.0)]
    for _ in range(k):
      lo = min(l for l, r in T)
      hi = max(r for l, r in T)
      delta = hi - lo
      S = []
      for start in (2, 6, 10, 14):
        offset = delta * start - lo
        for (l, r) in T:
          S.append((l + offset, r + offset))
      S += [
        (delta * 1,  delta * 5),
        (delta * 12, delta * 16),
        (delta * 4,  delta * 9),
        (delta * 8,  delta * 13),
      ]
      T = S
    return _normalize_grid(T)

  best_raw = best[5]
  best_omega = best[1]
  best_colors = best[2]
  best_ratio = best_colors / best_omega if best_omega > 0 else 0.0

  # Apply deterministic pruning to shrink the witness while preserving the frontier
  pruned = prune_by_ratio(best_raw, best_ratio)
  final = _normalize_grid(pruned)
  if not final:
    # fallback to normalized best (shouldn't happen, but safe)
    return _normalize_grid(best_raw)

  return final
# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()