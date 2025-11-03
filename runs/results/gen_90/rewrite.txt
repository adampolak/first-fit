# EVOLVE-BLOCK-START

def construct_intervals(iterations=4):
  """
  Construct a sequence of open intervals presented in order to FirstFit.
  Crossover implementation: flexible recursive builder + evaluation + pruning.

  Arguments:
    iterations: number of recursive expansion steps (default 4)

  Returns:
    normalized list of integer interval tuples (l, r)
  """

  # ------------------------------
  # Geometry & coloring utilities
  # ------------------------------
  def overlaps(a, b):
    """Open-interval overlap test: True iff intervals overlap."""
    (l1, r1), (l2, r2) = a, b
    return max(l1, l2) < min(r1, r2)

  def clique_number(intervals):
    """
    Compute omega (maximum number of intervals covering a single point) using sweep.
    Open intervals: process right(-1) before left(+1) at ties.
    """
    events = []
    for (l, r) in intervals:
      if l < r:
        events.append((l, +1))
        events.append((r, -1))
    if not events:
      return 0
    # sort by position, and close (-1) before open (+1) at same coordinate
    events.sort(key=lambda e: (e[0], 0 if e[1] == -1 else 1))
    cur = 0
    best = 0
    for _, t in events:
      cur += t
      if cur > best:
        best = cur
    return best

  def firstfit_colors(intervals):
    """
    Exact FirstFit simulator: iterate color classes and check overlaps.
    Correct for arbitrary arrival orders (not relying on start/end monotonicity).
    """
    colors = []  # list of lists of intervals assigned to each color (in arrival order)
    for iv in intervals:
      placed = False
      for c in colors:
        conflict = False
        for u in c:
          if overlaps(u, iv):
            conflict = True
            break
        if not conflict:
          c.append(iv)
          placed = True
          break
      if not placed:
        colors.append([iv])
    return len(colors)

  # ------------------------------
  # Normalization
  # ------------------------------
  def normalize_intervals(intervals):
    """
    Map unique endpoints to even integers starting from 0 (keeps order, avoids degenerate lengths).
    Returns list of integer (l, r) tuples.
    """
    if not intervals:
      return []
    points = sorted(set(x for seg in intervals for x in seg))
    coord = {}
    cur = 0
    for p in points:
      coord[p] = cur
      cur += 2
    return [(coord[l], coord[r]) for (l, r) in intervals]

  # ------------------------------
  # Flexible recursive pattern builder
  # ------------------------------
  def build_pattern(k, base_seed, offsets, blockers, translation='left',
                    blocker_anchor='left', schedule='after', interleave='block',
                    reverse_alt=False):
    """
    Build a k-level pattern:
      - base_seed: list of intervals for level 0
      - offsets: tuple of translation multipliers for copies at each level
      - blockers: list/tuple of (a,b) blocking intervals (in scaled coordinates)
      - translation: 'left' or 'center' anchor for copies
      - blocker_anchor: 'left' or 'center' anchor for blockers
      - schedule: 'after' | 'before' | 'split' decides ordering of copies vs blockers
      - interleave: 'block' (copies contiguously) or 'zip' (round-robin)
      - reverse_alt: if True, reverse order of every other copy to mix arrival order
    """
    T = list(base_seed)
    for _ in range(k):
      lo = min(l for l, r in T)
      hi = max(r for l, r in T)
      delta = hi - lo
      center = (lo + hi) / 2.0

      # create per-copy sequences
      copy_lists = []
      for idx, start in enumerate(offsets):
        if translation == 'left':
          off = delta * start - lo
        else:
          off = delta * start - center
        seq = T if not (reverse_alt and (idx % 2 == 1)) else list(reversed(T))
        copy_lists.append([(l + off, r + off) for (l, r) in seq])

      # interleave copies if requested
      if interleave == 'zip' and copy_lists:
        m = len(copy_lists[0])
        S_copies = []
        for j in range(m):
          for lst in copy_lists:
            S_copies.append(lst[j])
      else:
        S_copies = []
        for lst in copy_lists:
          S_copies.extend(lst)

      # scale blockers for this level
      S_blockers = []
      for (a, b) in blockers:
        if blocker_anchor == 'left':
          S_blockers.append((delta * a, delta * b))
        else:
          S_blockers.append((delta * a - center, delta * b - center))

      # compose arrival order according to schedule
      if schedule == 'before':
        S = S_blockers + S_copies
      elif schedule == 'split' and copy_lists:
        h = max(1, len(copy_lists) // 2)
        first_half = []
        for i in range(h):
          first_half.extend(copy_lists[i])
        second_half = []
        for i in range(h, len(copy_lists)):
          second_half.extend(copy_lists[i])
        S = first_half + S_blockers + second_half
      else:  # default 'after'
        S = S_copies + S_blockers

      T = S
    return T

  # ------------------------------
  # Evaluation + pruning helpers
  # ------------------------------
  def evaluate(raw_intervals):
    """
    Normalize raw_intervals and compute (score, omega, cols, n, normalized_intervals).
    Score is ratio with a tiny penalty for large n to prefer concise witnesses.
    """
    norm = normalize_intervals(raw_intervals)
    n = len(norm)
    if n == 0:
      return (-1.0, 0, 0, n, norm)
    om = clique_number(norm)
    if om == 0:
      return (-1.0, 0, 0, n, norm)
    cols = firstfit_colors(norm)
    ratio = cols / om
    score = ratio - 1e-6 * (n / 10000.0)
    return (score, om, cols, n, norm)

  def prune_strict(intervals, target_cols, target_omega):
    """
    Strict pruning: remove intervals (preferring longest) whenever both FF colors
    and clique number remain exactly equal to the target witness.
    """
    cur = list(intervals)
    def length(iv): return iv[1] - iv[0]
    order = sorted(range(len(cur)), key=lambda i: (-length(cur[i]), i))
    changed = True
    while changed:
      changed = False
      for idx in list(order):
        if idx >= len(cur):
          continue
        cand = cur[:idx] + cur[idx+1:]
        _, om, cols, _, _ = evaluate(cand)
        if om == target_omega and cols == target_cols:
          cur = cand
          order = sorted(range(len(cur)), key=lambda i: (-length(cur[i]), i))
          changed = True
          break
    return cur

  def prune_relaxed(intervals, min_ratio, max_omega):
    """
    Relaxed pruning: remove an interval if resulting ratio >= min_ratio
    and omega <= max_omega. Prefer to remove longer intervals first.
    """
    cur = list(intervals)
    def length(iv): return iv[1] - iv[0]
    order = sorted(range(len(cur)), key=lambda i: (-length(cur[i]), i))
    changed = True
    while changed:
      changed = False
      for idx in list(order):
        if idx >= len(cur):
          continue
        cand = cur[:idx] + cur[idx+1:]
        _, om, cols, _, _ = evaluate(cand)
        if om == 0 or om > max_omega:
          continue
        ratio = cols / om
        if ratio + 1e-12 >= min_ratio and om <= max_omega:
          cur = cand
          order = sorted(range(len(cur)), key=lambda i: (-length(cur[i]), i))
          changed = True
          break
    return cur

  # ------------------------------
  # Candidate generation (compact search)
  # ------------------------------
  offsets_set = [
    (2, 6, 10, 14),  # baseline
    (1, 5, 9, 13),   # shifted
    (3, 7, 11, 15),  # alternate
    (2, 6, 9, 14)    # skewed
  ]
  blockers_templates = [
    ((1, 5), (12, 16), (4, 9), (8, 13)),
    ((0, 4), (6, 10), (8, 12), (14, 18)),
    ((1, 6), (11, 16), (4, 8), (10, 14)),
  ]
  translations = ['left', 'center']
  blocker_anchors = ['left', 'center']
  depths = [max(2, iterations - 1), iterations]  # try a couple of depths
  schedules = ['after', 'split']
  interleaves = ['block', 'zip']
  reverse_flags = [False, True]
  extras = [False, True]  # optional extra copy on the first level

  best = None  # store tuple (score, om, cols, n, normalized, raw)
  # enumerate compact configuration space with a conservative size guard
  for depth in depths:
    for offsets in offsets_set:
      for blockers in blockers_templates:
        for translation in translations:
          for anchor in blocker_anchors:
            for schedule in schedules:
              for interleave in interleaves:
                for rev in reverse_flags:
                  for extra in extras:
                    # compute approximate size and skip exploding configs
                    offs = list(offsets)
                    if extra:
                      offs = offs + [max(offs) + 4]
                    copies_per_level = len(offs)
                    # rough final size: copies_per_level^depth (dominant) + blockers*depth
                    approx_size = (copies_per_level ** depth) * 1 + depth * len(blockers)
                    if approx_size > 4000:
                      continue

                    raw = build_pattern(
                      k=depth,
                      base_seed=[(0.0, 1.0)],
                      offsets=tuple(offs),
                      blockers=blockers,
                      translation=translation,
                      blocker_anchor=anchor,
                      schedule=schedule,
                      interleave=interleave,
                      reverse_alt=rev
                    )
                    score, om, cols, n, norm = evaluate(raw)
                    cand = (score, om, cols, n, norm, raw)
                    if best is None:
                      best = cand
                      continue
                    # pick best by score, tie-break by fewer intervals, then more colors
                    if cand[0] > best[0] + 1e-12:
                      best = cand
                    elif abs(cand[0] - best[0]) <= 1e-12:
                      if cand[3] < best[3]:
                        best = cand
                      elif cand[3] == best[3] and cand[2] > best[2]:
                        best = cand

  # Fallback canonical pattern (if search failed) - emulate Figure-4 baseline
  if best is None:
    T = [(0.0, 1.0)]
    for _ in range(4):
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
        (delta * 8,  delta * 13)
      ]
      T = S
    return normalize_intervals(T)

  # ------------------------------
  # Prune the best candidate (strict then relaxed), then final greedy shrink
  # ------------------------------
  _, best_om, best_cols, _, best_norm, best_raw = best

  # Strict pruning: keep exact witness
  phase1 = prune_strict(best_raw, best_cols, best_om)

  # Relaxed pruning: preserve ratio and omega cap
  target_ratio = best_cols / best_om if best_om > 0 else 0.0
  phase2 = prune_relaxed(phase1, target_ratio, best_om)

  # Final greedy shrink: try removing any interval that preserves ratio >= target_ratio and omega <= best_om
  final_raw = list(phase2)
  improved = True
  while improved:
    improved = False
    # try removing longer intervals first (heuristic)
    order = sorted(range(len(final_raw)), key=lambda i: -(final_raw[i][1] - final_raw[i][0]))
    for idx in order:
      cand = final_raw[:idx] + final_raw[idx+1:]
      score_c, om_c, cols_c, _, _ = evaluate(cand)
      if om_c == 0:
        continue
      if om_c <= best_om and cols_c / om_c + 1e-12 >= target_ratio:
        final_raw = cand
        improved = True
        break

  final_norm = normalize_intervals(final_raw)
  # safety: if ratio decreased unexpectedly, return the pre-pruned best normalized
  om_final = clique_number(final_norm)
  cols_final = firstfit_colors(final_norm)
  if om_final == 0 or cols_final / om_final + 1e-12 < target_ratio or om_final > best_om:
    return best_norm
  return final_norm

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()