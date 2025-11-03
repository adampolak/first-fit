# EVOLVE-BLOCK-START

def construct_intervals(iterations=4):
  """
  Construct a sequence of intervals of real line to maximize FirstFit/omega.

  Enhancements made here:
   - support interleave='zip' (round-robin emission of copies) and reverse_alt
     to increase FirstFit pressure without increasing clique size;
   - memoize FirstFit and clique computations to speed enumerations;
   - deterministic pruning that removes intervals while preserving both
     FirstFit colors and clique number (reduces instance size deterministically);
   - normalize once at the end of candidate construction (preserves geometry
     during recursion).
  """

  # --- Lightweight memoized metrics ---

  _ff_cache = {}
  _om_cache = {}

  def firstfit_colors(intervals):
    """Memoized FirstFit (arrival order given)."""
    key = tuple(intervals)
    if key in _ff_cache:
      return _ff_cache[key]
    last_end = []
    for (l, r) in intervals:
      placed = False
      for i, le in enumerate(last_end):
        # open intervals: if l >= le, can reuse color
        if l >= le:
          last_end[i] = r
          placed = True
          break
      if not placed:
        last_end.append(r)
    res = len(last_end)
    _ff_cache[key] = res
    return res

  def clique_number(intervals):
    """Memoized sweep-line for open intervals; right endpoints before left at ties."""
    key = tuple(intervals)
    if key in _om_cache:
      return _om_cache[key]
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
    _om_cache[key] = best
    return best

  def _normalize_grid(intervals):
    """Map unique endpoints to even integers (0,2,4,...)."""
    if not intervals:
      return []
    endpoints = sorted(set(x for seg in intervals for x in seg))
    coord = {}
    cur = 0
    for e in endpoints:
      coord[e] = cur
      cur += 2
    return [(coord[l], coord[r]) for (l, r) in intervals]

  # deterministic pruning: remove single intervals if both metrics unchanged
  def prune_preserve_metrics(intervals):
    if not intervals:
      return intervals
    # Work deterministically left-to-right, repeat until no change
    final = list(intervals)
    try:
      om = clique_number(final)
      cols = firstfit_colors(final)
      changed = True
      while changed:
        changed = False
        # iterate left-to-right deterministic indices
        i = 0
        while i < len(final):
          cand = final[:i] + final[i+1:]
          if not cand:
            i += 1
            continue
          if clique_number(cand) == om and firstfit_colors(cand) == cols:
            final = cand
            changed = True
            # do not increment i (we now examine the new element at i)
          else:
            i += 1
    except Exception:
      # if any error, return original
      return intervals
    return final

  def build_candidate(k, schedule='after', extra_first=False, translation='left',
                      offsets=(2,6,10,14), interleave='block', reverse_alt=False):
    """
    Build k-level recursive pattern using copies and blockers.
    interleave in {'block', 'zip'}:
      - 'block' emits all intervals of the first copy, then second, ...
      - 'zip' emits one interval from each copy in round-robin order (zip).
    reverse_alt: if True and interleave='block', every other copy's internal order
      is reversed before emission; for 'zip' it reverses the copy sequence to break symmetry.
    """
    T = [(0.0, 1.0)]
    for i in range(k):
      lo = min(l for l, r in T)
      hi = max(r for l, r in T)
      delta = hi - lo
      center = (lo + hi) / 2.0

      # optionally include a fifth copy on the first iteration
      offs = tuple(list(offsets) + ([18] if extra_first and i == 0 else []))

      # Build per-copy lists (without flattening)
      copy_lists = []
      for idx, start in enumerate(offs):
        if translation == 'left':
          offset = delta * start - lo
        else:
          offset = delta * start - center
        # optionally reverse internal order for alternating copies to break alignment
        seq = list(reversed(T)) if (reverse_alt and (idx % 2 == 1)) else list(T)
        copy_lists.append([(l + offset, r + offset) for (l, r) in seq])

      # Interleave options
      if interleave == 'block':
        # flatten copies in block order
        S_copies = []
        for lst in copy_lists:
          S_copies.extend(lst)
      else:  # 'zip' - round-robin emission
        S_copies = []
        if copy_lists:
          m = max(len(lst) for lst in copy_lists)
          for j in range(m):
            for lst in copy_lists:
              if j < len(lst):
                S_copies.append(lst[j])

      # Blockers (left-anchored to maintain canonical geometry)
      blockers = [
        (delta * 1,  delta * 5),
        (delta * 12, delta * 16),
        (delta * 4,  delta * 9),
        (delta * 8,  delta * 13),
      ]

      if schedule == 'after':
        S = S_copies + blockers
      elif schedule == 'before':
        S = blockers + S_copies
      else:  # split: half copies, blockers, remainder
        h = max(1, len(copy_lists) // 2)
        first = []
        for lst in copy_lists[:h]:
          first.extend(lst)
        second = []
        for lst in copy_lists[h:]:
          second.extend(lst)
        S = first + blockers + second

      # keep full-precision geometry during recursion (normalize later)
      T = S

    # optionally run a light deterministic pruning (keeps metrics unchanged)
    T_pruned = prune_preserve_metrics(T)

    # normalize once before returning/evaluation
    return _normalize_grid(T_pruned)

  # --- Enumerate a modest but richer set of candidates and select the best ---

  depths = [max(2, iterations - 1), iterations]
  schedules = ['after', 'split']  # 'before' rarely helps; keep search tight
  extras = [False, True]
  translations = ['left', 'center']
  interleaves = ['block', 'zip']
  rev_flags = [False, True]
  offsets_options = [
    (2, 6, 10, 14),      # baseline
    (1, 5, 9, 13),       # shifted
    (3, 7, 11, 15),      # alternate
    (2, 6, 10, 14, 18),  # extra copy variant
    (1.5, 5.5, 9.5, 13.5) # half-step stagger
  ]

  best_T = None
  best_ratio = -1.0
  best_n = None
  best_cols = 0
  best_om = 0

  for k in depths:
    for sch in schedules:
      for ex in extras:
        for trans in translations:
          for inter in interleaves:
            for rev in rev_flags:
              for offs in offsets_options:
                # build candidate
                T = build_candidate(k, schedule=sch, extra_first=ex,
                                    translation=trans, offsets=offs,
                                    interleave=inter, reverse_alt=rev)
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

  # Fallback to canonical baseline if none found (shouldn't happen)
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
    return _normalize_grid(T)

  # Final pruning on the chosen best (deterministic, metric-preserving)
  final_T = prune_preserve_metrics(best_T)
  return _normalize_grid(final_T)

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()