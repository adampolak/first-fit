# EVOLVE-BLOCK-START

from typing import List, Tuple

Interval = Tuple[float, float]

def construct_intervals(iterations=4) -> List[Tuple[int, int]]:
  """
  Construct a sequence of open intervals of the real line, in arrival order,
  to maximize the FirstFit/omega ratio. We blend the canonical Figure-4
  recursive skeleton with:
    - rotor tilings across levels (four offset sets rotated deterministically),
    - alternate blocker templates,
    - deterministic, frontier-preserving pruning,
    - endpoint normalization to a compact even-integer grid.

  Input:
    iterations: recursion depth (default 4).

  Output:
    list of (l, r) integer endpoints with l < r (open intervals).
  """

  # ---------- Utilities ----------
  def overlaps(a: Interval, b: Interval) -> bool:
    # Open intervals overlap iff not (a.r <= b.l or b.r <= a.l)
    (al, ar), (bl, br) = a, b
    return not (ar <= bl or br <= al)

  def firstfit_colors_and_assignment(intervals: List[Interval]):
    # Full FirstFit: smallest color whose class has no overlapping interval.
    color_classes: List[List[Interval]] = []
    assignment: List[int] = []
    for seg in intervals:
      placed = False
      for c, cls in enumerate(color_classes):
        if all(not overlaps(seg, x) for x in cls):
          cls.append(seg)
          assignment.append(c)
          placed = True
          break
      if not placed:
        color_classes.append([seg])
        assignment.append(len(color_classes) - 1)
    return len(color_classes), assignment

  def clique_number(intervals: List[Interval]) -> int:
    # Exact sweep for open intervals: at ties, process -1 before +1.
    events = []
    for l, r in intervals:
      if l < r:
        events.append((l, +1))
        events.append((r, -1))
    events.sort(key=lambda e: (e[0], e[1]))  # (-1) before (+1) at ties
    cur = best = 0
    for _, delta in events:
      cur += delta
      if cur > best:
        best = cur
    return best

  def normalize_even_grid(intervals: List[Interval]) -> List[Tuple[int, int]]:
    if not intervals:
      return []
    pts = sorted({x for seg in intervals for x in seg})
    mapping = {x: 2*i for i, x in enumerate(pts)}
    return [(mapping[l], mapping[r]) for l, r in intervals]

  # ---------- Geometry blueprint ----------
  starts_rotor = [
    (2, 6, 10, 14),  # A
    (1, 5, 9, 13),   # B
    (3, 7, 11, 15),  # C
    (0, 4, 8, 12),   # D
  ]
  blocker_sets = [
    [(1, 5), (12, 16), (4, 9), (8, 13)],  # baseline
    [(1, 6), (11, 16), (3, 9), (7, 13)],  # elongated/shifted
    [(2, 6), (12, 16), (4, 9), (8, 13)],  # shifted 1st
  ]

  def expand_once(T: List[Interval], offsets: Tuple[float, ...], blockers_ab: List[Tuple[float, float]]) -> List[Interval]:
    lo = min(l for l, _ in T)
    hi = max(r for _, r in T)
    d = hi - lo
    S: List[Interval] = []
    for o in offsets:
      off = d * o - lo
      for l, r in T:
        S.append((l + off, r + off))
    for a, b in blockers_ab:
      S.append((d * a, d * b))
    return S

  def generate(depth: int, use_rotor: bool, blocker_template_idx: int, extra_first_copy: bool) -> List[Interval]:
    T: List[Interval] = [(0.0, 1.0)]
    for lvl in range(depth):
      if use_rotor:
        offsets = starts_rotor[lvl % 4]
      else:
        offsets = starts_rotor[0]
      # optional deterministic augmentation: add a 5th copy on the first level
      if extra_first_copy and lvl == 0:
        offsets = tuple(list(offsets) + [18])
      S = expand_once(T, offsets, blocker_sets[blocker_template_idx])
      T = S
    return T

  # ---------- Deterministic pruning preserving ratio ----------
  def prune_preserve_ratio(intervals: List[Interval], keep_ratio: float, max_rounds: int = 2, budget_per_round: int = 200) -> List[Interval]:
    T = list(intervals)
    # Stage order: first remove longer intervals, then general pass.
    for round_id in range(max_rounds):
      attempted = 0
      changed = True
      while changed and attempted < budget_per_round:
        changed = False
        if round_id == 0:
          order = sorted(range(len(T)), key=lambda i: -(T[i][1] - T[i][0]))
        else:
          order = list(range(len(T)))
        for idx in order:
          if attempted >= budget_per_round:
            break
          attempted += 1
          cand = T[:idx] + T[idx+1:]
          om = clique_number(cand)
          if om == 0:
            continue
          cols, _ = firstfit_colors_and_assignment(cand)
          if cols / om >= keep_ratio - 1e-12:
            T = cand
            changed = True
            break
    return T

  # ---------- Candidate search ----------
  # Ensure we always include the canonical baseline (depth=iterations, rotor off, blockers idx 0)
  candidates = []
  param_grid = []
  depths = [max(2, iterations - 1), iterations]
  for depth in depths:
    for rotor_flag in (False, True):
      for blk_idx in range(len(blocker_sets)):
        for extra in (False, True):
          param_grid.append((depth, rotor_flag, blk_idx, extra))

  best_T = None
  best_ratio = -1.0
  best_len = None
  for depth, rotor_flag, blk_idx, extra in param_grid:
    raw = generate(depth, rotor_flag, blk_idx, extra)
    # Evaluate pre-normalization (affine transforms; normalization preserves overlaps)
    om = clique_number(raw)
    if om == 0:
      continue
    cols, _ = firstfit_colors_and_assignment(raw)
    ratio = cols / om
    # Keep the best; tie-breaker: shorter instance
    if (ratio > best_ratio + 1e-12) or (abs(ratio - best_ratio) <= 1e-12 and (best_len is None or len(raw) < best_len)):
      best_T = raw
      best_ratio = ratio
      best_len = len(raw)

  # Fallback: canonical if search failed
  if best_T is None:
    best_T = generate(iterations, False, 0, False)
    om = clique_number(best_T)
    cols, _ = firstfit_colors_and_assignment(best_T)
    best_ratio = cols / max(1, om)

  # Conservative pruning: keep the current ratio
  pruned = prune_preserve_ratio(best_T, best_ratio, max_rounds=2, budget_per_round=100)

  # Final normalization
  return normalize_even_grid(pruned)

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()