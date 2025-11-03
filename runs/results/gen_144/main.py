# EVOLVE-BLOCK-START

from typing import List, Tuple

Interval = Tuple[int, int]

def _normalize_grid(intervals: List[Tuple[float, float]]) -> List[Interval]:
  """
  Normalize endpoints to a compact integer grid while preserving order.
  Each unique endpoint is mapped to an increasing even integer.
  Returns a new list of (l, r) with integer coordinates.
  """
  if not intervals:
    return []
  endpoints = sorted(set([x for seg in intervals for x in seg]))
  coord = {}
  cur = 0
  for e in endpoints:
    coord[e] = cur
    cur += 2  # spacing by 2 to keep even gaps

  return [(coord[l], coord[r]) for (l, r) in intervals]


def _overlaps(iv1: Interval, iv2: Interval) -> bool:
  """Open-interval overlap test: return True iff intervals overlap."""
  (l1, r1), (l2, r2) = iv1, iv2
  return max(l1, l2) < min(r1, r2)


def _clique_number(intervals: List[Interval]) -> int:
  """
  Compute omega (maximum number of intervals covering a single point) using sweep.
  For open intervals, at identical coordinates, process right endpoints before left.
  """
  events = []
  for (l, r) in intervals:
    if l < r:
      events.append((l, +1))
      events.append((r, -1))
  # Right endpoints first at same coordinate for open intervals
  events.sort(key=lambda e: (e[0], 0 if e[1] == -1 else 1))
  cur = 0
  best = 0
  for _, t in events:
    cur += t
    if cur > best:
      best = cur
  return best


def _firstfit_colors(intervals: List[Interval]) -> int:
  """
  Simulate FirstFit coloring on the given arrival order.
  Return total number of colors used. O(n^2) exact simulation.
  """
  colors: List[List[Interval]] = []
  for iv in intervals:
    placed = False
    for c in colors:
      # Color class invariant: pairwise non-overlapping; check against last may suffice,
      # but we check all for robustness with arbitrary orders.
      conflict = False
      for u in c:
        if _overlaps(u, iv):
          conflict = True
          break
      if not conflict:
        c.append(iv)
        placed = True
        break
    if not placed:
      colors.append([iv])
  return len(colors)


def _make_copies(T, offsets, delta, lo, center, translation):
  """
  Create translated copies of T according to offsets and translation rule.
  translation in {'left', 'center'}.
  """
  S = []
  for start in offsets:
    if translation == 'left':
      offset = delta * start - lo
    else:  # center-based
      offset = delta * start - center
    for (l, r) in T:
      S.append((l + offset, r + offset))
  return S


def _add_blockers(S, blockers, delta, anchor, center):
  """
  Add blockers scaled by delta. Anchor may be 'left' (absolute) or 'center' (center-shifted).
  """
  for (a, b) in blockers:
    if anchor == 'left':
      S.append((delta * a, delta * b))
    else:
      S.append((delta * a - center, delta * b - center))
  return S


def _build_pattern(iterations, base_seed, tiling_schedule, blockers, translation, blocker_anchor):
  """
  Recursively expand the base_seed 'iterations' times using a possibly rotating
  tiling schedule (sequence of offset tuples), adding blockers at each level.
  """
  T = list(base_seed)
  for i in range(iterations):
    lo = min(l for l, r in T)
    hi = max(r for l, r in T)
    delta = hi - lo
    center = (lo + hi) / 2.0

    offsets = tiling_schedule[i % len(tiling_schedule)]
    S = _make_copies(T, offsets, delta, lo, center, translation)
    S = _add_blockers(S, blockers, delta, blocker_anchor, center)
    T = S
  return T


def _prune_preserve(intervals: List[Interval]) -> List[Interval]:
  """
  Deterministic pruning: iteratively remove intervals if both
  - omega (clique number) stays the same, and
  - FirstFit color count stays the same.
  This tends to reduce footprint with the same ratio.
  """
  if not intervals:
    return intervals
  base_cols = _firstfit_colors(intervals)
  base_omega = _clique_number(intervals)
  n = len(intervals)
  # Try removals in a deterministic order: by descending length, then by left, then by right
  order = sorted(range(n), key=lambda i: (-(intervals[i][1] - intervals[i][0]), intervals[i][0], intervals[i][1]))
  kept = list(intervals)
  removed_any = True
  while removed_any:
    removed_any = False
    # Recompute order for current kept
    order = sorted(range(len(kept)), key=lambda i: (-(kept[i][1] - kept[i][0]), kept[i][0], kept[i][1]))
    for idx in order:
      cand = kept[:idx] + kept[idx+1:]
      if not cand:
        continue
      om = _clique_number(cand)
      cols = _firstfit_colors(cand)
      if om == base_omega and cols == base_cols:
        kept = cand
        removed_any = True
        break
  return kept


def construct_intervals(iterations=4, normalize=True):
  """
  Build a sequence of open intervals that forces FirstFit
  to use many colors while keeping the clique size controlled.
  This generator explores several deterministic blueprints based on
  the 4-copy + 4-blocker scheme and selects the best candidate.
  """
  # Base seeds to try
  base_seeds = [
    [(0.0, 1.0)],  # canonical
  ]

  # Offset tilings (four copies per level), rotated as a schedule
  A = (2, 6, 10, 14)
  B = (1, 5, 9, 13)
  C = (3, 7, 11, 15)
  D = (0, 4, 8, 12)

  tiling_schedules = [
    (A,),           # fixed A
    (B,),           # fixed B
    (C,),           # fixed C
    (D,),           # fixed D
    (A, B, C, D),   # rotating A->B->C->D
    (A, C, B, D),   # alternate rotation
    (A, D, B, C),   # another rotation
  ]

  # Blocker templates
  blockers_templates = [
    ((1, 5), (12, 16), (4, 9), (8, 13)),   # baseline A
    ((0, 4), (6, 10), (8, 12), (14, 18)),  # template B
    ((2, 6), (4, 8), (10, 14), (12, 16)),  # template C
  ]

  translations = ['left', 'center']
  anchors = ['left', 'center']

  # Evaluate all candidates and choose the best per ratio; tie-break by smaller n.
  best_norm: List[Interval] = []
  best_ratio = float('-inf')
  best_cols = -1
  best_omega = 1
  best_n = 10**9

  for seed in base_seeds:
    for schedule in tiling_schedules:
      for blockers in blockers_templates:
        for translation in translations:
          for anchor in anchors:
            # Build raw geometry
            raw = _build_pattern(iterations, seed, schedule, blockers, translation, anchor)
            # Normalize to integer grid if requested
            Tn = _normalize_grid(raw) if normalize else [(l, r) for (l, r) in raw]
            if not Tn:
              continue
            om = _clique_number(Tn)
            if om <= 0:
              continue
            cols = _firstfit_colors(Tn)
            ratio = cols / om
            nlen = len(Tn)

            # Deterministic pruning stage to reduce footprint without hurting metrics
            pruned = _prune_preserve(Tn)
            om_p = _clique_number(pruned)
            cols_p = _firstfit_colors(pruned)
            ratio_p = cols_p / om_p if om_p > 0 else float('-inf')
            # Accept pruning if it keeps the ratio and color/omega
            if ratio_p >= ratio and cols_p >= cols and om_p == om:
              Tn = pruned
              ratio = ratio_p
              cols = cols_p
              nlen = len(Tn)

            # Keep best by ratio; tie-break by fewer intervals, then more colors
            if ratio > best_ratio + 1e-12 or (abs(ratio - best_ratio) <= 1e-12 and (nlen < best_n or (nlen == best_n and cols > best_cols))):
              best_ratio = ratio
              best_cols = cols
              best_omega = om
              best_n = nlen
              best_norm = Tn

  # Fallback: if search produced nothing, return classic construction
  if not best_norm:
    # classic four-copy construction
    T = [(0.0, 1.0)]
    for _ in range(iterations):
      lo = min(l for l, r in T)
      hi = max(r for l, r in T)
      delta = hi - lo
      S = []
      for start in (2, 6, 10, 14):
        offset = delta * start - lo
        S.extend([(offset + l, offset + r) for l, r in T])
      S += [
        (delta * 1,  delta * 5),
        (delta * 12, delta * 16),
        (delta * 4,  delta * 9),
        (delta * 8,  delta * 13)
      ]
      T = S
    return _normalize_grid(T) if normalize else T

  return best_norm

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()