# EVOLVE-BLOCK-START

def construct_intervals(iterations: int = 4):
  """
  Construct a sequence of intervals on the real line, presented in order to FirstFit,
  aiming to maximize the ratio FirstFit_colors / omega (max clique size).
  We perform a compact search over small adversarial variations:
    - translation: 'left' vs 'center'
    - offsets pattern: several stagger placements
    - schedule of blockers (before/after/split)
    - optional extra copy on the first level
    - weaving copies (round-robin) to couple colors across substructures
    - optional reversal at odd levels (arrival-order engineering)
  Then select the candidate maximizing the simulated ratio and normalize to a
  compact even-integer grid, preserving intersection relations for open intervals.

  Args:
    iterations: recursion depth (default 4)

  Returns:
    List[Tuple[int,int]]: open intervals (l, r) with l < r, on an even integer grid.
  """

  # ---- Helpers ----

  def firstfit_colors(intervals):
    """
    Simulate FirstFit on open intervals. For each color, keep the last end; reuse if
    new interval starts >= that end (open interval model).
    """
    last_end = []
    for (l, r) in intervals:
      placed = False
      for i, e in enumerate(last_end):
        if l >= e:
          last_end[i] = r
          placed = True
          break
      if not placed:
        last_end.append(r)
    return len(last_end)

  def clique_number(intervals):
    """
    Compute omega using sweep over open intervals: process end events before start
    at equal positions to model openness.
    """
    ev = []
    for (l, r) in intervals:
      if l < r:
        ev.append((l, +1))
        ev.append((r, -1))
    # At ties, end (-1) before start (+1)
    ev.sort(key=lambda x: (x[0], 0 if x[1] == -1 else 1))
    cur = best = 0
    for _, d in ev:
      cur += d
      if cur > best:
        best = cur
    return best

  def normalize_even_grid(intervals):
    """
    Order-preserving remap of all endpoints to even integers (0,2,4,...).
    """
    if not intervals:
      return []
    endpoints = sorted({x for seg in intervals for x in seg})
    m = {}
    cur = 0
    for e in endpoints:
      m[e] = cur
      cur += 2
    return [(m[l], m[r]) for (l, r) in intervals]

  def build_candidate(
      k,
      schedule='after',         # 'after', 'before', 'split'
      extra_first=False,        # add a 5th copy (offset 18) only on level 0
      translation='left',       # 'left' or 'center'
      offsets=(2, 6, 10, 14),
      weave=False,              # interleave copies round-robin
      reverse_odd_levels=False  # reverse arrival order on odd levels
  ):
    """
    Build a k-level recursive adversary.

    - T starts as [(0,1)].
    - Each level: create several translated copies plus four long 'blockers'.
    - Schedule controls arrival order of copies vs blockers.
    - Weave interleaves the copies to increase coupling across colors.
    - reverse_odd_levels optionally reverses order at odd levels (0-indexed).
    """
    T = [(0, 1)]
    cur_lo, cur_hi = 0, 1  # track current extent
    delta = cur_hi - cur_lo

    for lvl in range(k):
      # determine offsets for this level
      offs = list(offsets)
      if extra_first and lvl == 0:
        offs.append(18)

      # choose anchor for translation
      anchor = cur_lo if translation == 'left' else (cur_lo + cur_hi) / 2.0

      # Prepare the copies
      # Optionally weave them: round-robin one interval from each translated copy
      def make_copies(subset):
        nonlocal delta, anchor, T
        res = []
        min_l, max_r = float('inf'), -float('inf')
        if not subset:
          return res, min_l, max_r

        if weave:
          shifted_lists = []
          for s in subset:
            shift = delta * s - anchor
            shifted_lists.append([(l + shift, r + shift) for (l, r) in T])
          # round-robin weave
          for i in range(len(T)):
            for arr in shifted_lists:
              l, r = arr[i]
              res.append((l, r))
              if l < min_l: min_l = l
              if r > max_r: max_r = r
        else:
          for s in subset:
            shift = delta * s - anchor
            for (l, r) in T:
              L, R = (l + shift, r + shift)
              res.append((L, R))
              if L < min_l: min_l = L
              if R > max_r: max_r = R
        return res, min_l, max_r

      # Four long connectors (blockers)
      blockers_raw = [
        (delta * 1,  delta * 5),
        (delta * 12, delta * 16),
        (delta * 4,  delta * 9),
        (delta * 8,  delta * 13),
      ]

      # Build S and track extents during construction to avoid extra scans
      S = []
      minL, maxR = float('inf'), -float('inf')

      if schedule == 'after':
        copies, aL, aR = make_copies(offs)
        S.extend(copies)
        if aL < minL: minL = aL
        if aR > maxR: maxR = aR
        for seg in blockers_raw:
          S.append(seg)
          l, r = seg
          if l < minL: minL = l
          if r > maxR: maxR = r
      elif schedule == 'before':
        for seg in blockers_raw:
          S.append(seg)
          l, r = seg
          if l < minL: minL = l
          if r > maxR: maxR = r
        copies, aL, aR = make_copies(offs)
        S.extend(copies)
        if aL < minL: minL = aL
        if aR > maxR: maxR = aR
      else:  # 'split'
        h = max(1, len(offs) // 2)
        first, second = offs[:h], offs[h:]
        c1, l1, r1 = make_copies(first)
        S.extend(c1)
        if l1 < minL: minL = l1
        if r1 > maxR: maxR = r1
        for seg in blockers_raw:
          S.append(seg)
          l, r = seg
          if l < minL: minL = l
          if r > maxR: maxR = r
        c2, l2, r2 = make_copies(second)
        S.extend(c2)
        if l2 < minL: minL = l2
        if r2 > maxR: maxR = r2

      # Optional reversal at odd levels (arrival-order engineering)
      if reverse_odd_levels and (lvl % 2 == 1):
        S.reverse()

      # Next level
      T = S
      cur_lo, cur_hi = minL, maxR
      delta = cur_hi - cur_lo

    return T

  # ---- Small configuration search space ----

  depths = [max(2, iterations - 1), iterations]
  schedules = ['after', 'before', 'split']
  translations = ['left', 'center']
  offsets_options = [
    (2, 6, 10, 14),   # classic
    (1, 5, 9, 13),    # shifted
    (3, 7, 11, 15),   # alternate
    (0, 4, 8, 12),    # dense near-left
  ]
  extras = [False, True]
  weaves = [False, True]
  reverse_odd = [False, True]

  best_T = None
  best_ratio = -1.0
  best_n = None
  best_cols = 0
  best_om = 0

  for k in depths:
    for sch in schedules:
      for trans in translations:
        for offs in offsets_options:
          for ex in extras:
            for w in weaves:
              for rev in reverse_odd:
                T = build_candidate(
                  k,
                  schedule=sch,
                  extra_first=ex,
                  translation=trans,
                  offsets=offs,
                  weave=w,
                  reverse_odd_levels=rev
                )
                om = clique_number(T)
                if om == 0:
                  continue
                cols = firstfit_colors(T)
                ratio = cols / om
                n = len(T)

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

  # Fallback to baseline construction if search failed for some reason
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
    return normalize_even_grid(T)

  # Return the best candidate normalized to an even integer grid
  return normalize_even_grid(best_T)

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()