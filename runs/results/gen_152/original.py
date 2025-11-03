# EVOLVE-BLOCK-START

def construct_intervals(iterations=4):
  """
  Build a sequence of open intervals that aims to maximize FirstFit/OPT.
  Strategy:
    - 4-copy recursive expansion + 4 blockers per level (Figure-4 style)
    - Parameter sweep over offsets, blocker templates, translation anchors, and depth
    - Evaluate each candidate via an internal FirstFit simulation and clique sweep
    - Normalize endpoints to a compact integer grid before returning

  Arguments:
    iterations: optional hint; included in the depth sweep if in {3,4,5}

  Returns:
    intervals: list[(l,r)] of open intervals with integer endpoints
  """

  # ------------ helpers ------------
  def overlaps(a, b):
    (l1, r1), (l2, r2) = a, b
    # open-interval overlap:
    return max(l1, l2) < min(r1, r2)

  def firstfit_colors(intervals):
    """
    Simulate FirstFit on the given arrival order.
    Colors are independent sets w.r.t. open-interval overlap.
    """
    colors = []  # list of color classes
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

  def clique_number(intervals):
    """
    Max number of intervals covering a single point (omega) using sweep.
    For open intervals, at equal coordinates process right(-1) before left(+1).
    """
    events = []  # (x, type) type=-1 for right, +1 for left
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

  def normalize_grid(intervals):
    """
    Map each unique endpoint to an increasing integer grid (even spacing).
    Preserves open-interval overlap and arrival order.
    """
    if not intervals:
      return []
    pts = sorted(set([x for seg in intervals for x in seg]))
    coord = {e: 2*i for i, e in enumerate(pts)}
    return [(coord[l], coord[r]) for (l, r) in intervals]

  def build_pattern(base_seed, k, offsets, blockers, translation='left', blocker_anchor='left',
                    schedule='after', interleave='block', reverse_alt=False):
    """
    Recursively expand base_seed k times with multiple copies and blockers.
    - translation in {'left','center'} controls copy offsets
    - blocker_anchor in {'left','center'} controls blocker placement
    - schedule in {'after','before','split'} decides the order of copies vs blockers per level
    - interleave in {'block','zip'} controls how copies are interleaved
    - reverse_alt: if True and interleave='block', every other copy is presented in reverse order
    """
    T = list(base_seed)
    for _ in range(k):
      lo = min(l for l, r in T)
      hi = max(r for l, r in T)
      delta = hi - lo
      center = (lo + hi) / 2.0

      # Build per-copy sequences
      copy_lists = []
      for idx, start in enumerate(offsets):
        if translation == 'left':
          off = delta * start - lo
        else:  # center-based
          off = delta * start - center
        seq = T if not (reverse_alt and (idx % 2 == 1)) else list(reversed(T))
        copy_lists.append([(l + off, r + off) for (l, r) in seq])

      # Interleave copies
      if interleave == 'block':
        S_copies = []
        for lst in copy_lists:
          S_copies.extend(lst)
      else:  # 'zip': emit one item from each copy in round-robin order
        S_copies = []
        if copy_lists:
          m = len(copy_lists[0])
          for j in range(m):
            for lst in copy_lists:
              S_copies.append(lst[j])

      # Build blockers for this level
      S_blockers = []
      for (a, b) in blockers:
        if blocker_anchor == 'left':
          S_blockers.append((delta * a, delta * b))
        else:
          S_blockers.append((delta * a - center, delta * b - center))

      # Compose level by schedule
      if schedule == 'before':
        S = S_blockers + S_copies
      elif schedule == 'split' and interleave == 'block':
        h = max(1, len(copy_lists) // 2)
        first_half = []
        for i in range(h):
          first_half.extend(copy_lists[i])
        second_half = []
        for i in range(h, len(copy_lists)):
          second_half.extend(copy_lists[i])
        S = first_half + S_blockers + second_half
      elif schedule == 'mix':  # alternate copy then blocker intervals
        S = []
        max_len = max(len(S_copies), len(S_blockers))
        for j in range(max_len):
          if j < len(S_copies):
            S.append(S_copies[j])
          if j < len(S_blockers):
            S.append(S_blockers[j])
      else:
        # default 'after' or 'split' with interleave='zip'
        S = S_copies + S_blockers

      T = S
    return T

  def evaluate(intervals):
    """Return (score, omega, colors, n)."""
    if not intervals:
      return (-1.0, 0, 0, 0)
    om = clique_number(intervals)
    if om <= 0:
      return (-1.0, 0, 0, len(intervals))
    cols = firstfit_colors(intervals)
    # Slight penalty for larger instances to break ties
    score = cols / om - 1e-6 * (len(intervals) / 10000.0)
    return (score, om, cols, len(intervals))

  # ------------ parameter sweep ------------
  offsets_set = [
    (2, 6, 10, 14),      # canonical
    (1, 5, 9, 13),
    (3, 7, 11, 15),
    (0, 4, 8, 12),
    (2, 6, 9, 14),       # slight skew to tighten mid overlap
    (2, 7, 11, 14),      # shifted second copy
    (1.5, 5.5, 9.5, 13.5),  # half-step stagger
    (2, 6, 10, 14, 18)   # extra-copy offsets (explore 5-copy variants)
  ]
  blockers_templates = [
    # A: baseline as in Figure 4
    ((1, 5), (12, 16), (4, 9), (8, 13)),
    # B: variant geometry
    ((0, 4), (6, 10), (8, 12), (14, 18)),
    # C: tighter mid coupling
    ((2, 6), (4, 8), (10, 14), (12, 16)),
    # D: asymmetric early/late pushers
    ((1, 5), (12, 16), (3, 7), (9, 13)),
    # E: lengthened extremes
    ((1, 6), (11, 16), (4, 8), (10, 14)),
    # F: six-blocker couplers (two extra caps)
    ((1, 5), (12, 16), (4, 9), (8, 13), (2, 14), (6, 10)),
    # G: dense 5-blocker coupler (extra cap to more tightly couple copies)
    ((1, 5), (3, 9), (6, 10), (11, 15), (8, 13))
  ]
  translations = ['left', 'center']
  blocker_anchors = ['left', 'center']
  depths = [d for d in {3, 4, 5, iterations} if d in (3, 4, 5)]
  base_seeds = [
    [(0.0, 1.0)],
    [(0.0, 1.0), (2.0, 3.0)],
    [(0.0, 1.0), (0.5, 1.5)],  # new overlapping base seed to increase initial FF pressure
  ]

  best = None  # (score, om, cols, n, intervals)

  schedules = ['after', 'before', 'split', 'mix']
  interleaves = ['block', 'zip']
  rev_flags = [False, True]

  for base in base_seeds:
    for k in depths:
      # avoid overly large instances (approx size ~ 4^k * |base| + O(4^k))
      if (4 ** k) * (len(base) + 2) > 3000:
        continue
      for offsets in offsets_set:
        for blockers in blockers_templates:
          for tr in translations:
            for ba in blocker_anchors:
              for sch in schedules:
                for inter in interleaves:
                  for rev in rev_flags:
                    T = build_pattern(base_seed=base, k=k, offsets=offsets,
                                      blockers=blockers, translation=tr, blocker_anchor=ba,
                                      schedule=sch, interleave=inter, reverse_alt=rev)
                    score, om, cols, n = evaluate(T)
                    cand = (score, om, cols, n, T)
                    if best is None:
                      best = cand
                    else:
                      if cand[0] > best[0] + 1e-9:
                        best = cand
                      elif abs(cand[0] - best[0]) <= 1e-9:
                        # tie-break by fewer intervals, then more colors
                        if cand[3] < best[3] or (cand[3] == best[3] and cand[2] > best[2]):
                          best = cand

  # Fallback to the known-safe baseline if sweep failed (should not happen)
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
    return normalize_grid(T)

  # Attempt deterministic pruning: try removing single intervals (in arrival order)
  # if both the clique number and FirstFit color count remain unchanged.
  # This is conservative (single removals only) and preserves the adversarial effect.
  final_T = list(best[4])
  try:
    base_om = clique_number(final_T)
    base_cols = firstfit_colors(final_T)
    i = 0
    # single pass deterministic shrink: iterate once left-to-right
    while i < len(final_T):
      cand = final_T[:i] + final_T[i+1:]
      if not cand:
        i += 1
        continue
      om_c = clique_number(cand)
      cols_c = firstfit_colors(cand)
      # only accept removal if both metrics unchanged
      if om_c == base_om and cols_c == base_cols:
        # accept removal (do not increment i because list shifted left)
        final_T = cand
        # do not change base_om/base_cols
      else:
        i += 1
  except Exception:
    # If anything goes wrong in pruning, fall back to the unpruned best
    final_T = list(best[4])

  # Normalize endpoints to compact integer grid and return
  return normalize_grid(final_T)

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()