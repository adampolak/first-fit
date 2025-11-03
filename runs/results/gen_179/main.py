# EVOLVE-BLOCK-START

def construct_intervals(iterations=4):
  """
  Construct a sequence of intervals of real line to maximize FirstFit/omega.
  We explore small variations inspired by adversarial arrival-order engineering:
  - schedule of blockers vs copies at each level (before/after/split/interleaved)
  - optional extra copy at the first iteration to couple layers
  - center-anchored copy translations to break symmetry
  - compact integer normalization of endpoints before returning
  We evaluate candidates via an in-block FirstFit and omega computation, and return the best.
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
    events.sort(key=lambda e: (e[0], 0 if e[1] == -1 else 1))
    cur = best = 0
    for _, t in events:
      cur += t
      if cur > best:
        best = cur
    return best

  def _normalize_grid(intervals):
    """
    Map all unique endpoints to increasing even integers (0,2,4,...).
    This keeps instances compact and avoids unwieldy coordinates.
    Arrival order in the returned list is preserved.
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

  def build_candidate(k, schedule='after', extra_first=False, translation='left', offsets=(2, 6, 10, 14)):
    """
    Build k-level recursive pattern using scaled copies per level,
    with four long connectors ("blockers"), and a schedule controlling
    the arrival order within each level.

    translation in {'left','center'} controls whether copies/blockers
    are anchored by the current left endpoint or by the center.
    offsets may contain fractional values (e.g., half-step) to break symmetry.
    """
    T = [(0.0, 1.0)]
    for i in range(k):
      lo = min(l for l, r in T)
      hi = max(r for l, r in T)
      delta = hi - lo
      center = (lo + hi) / 2.0

      # anchor choice: left or center
      anchor = lo if translation == 'left' else center

      # optionally include a fifth copy on the first iteration
      offs = tuple(list(offsets) + ([18] if extra_first and i == 0 else []))

      def make_copies(from_T, offsets_local):
        S = []
        for start in offsets_local:
          offset = delta * start - anchor
          for (l, r) in from_T:
            S.append((l + offset, r + offset))
        return S

      # Blocker template (kept canonical); positioned relative to chosen anchor
      blocker_template = [(1, 5), (12, 16), (4, 9), (8, 13)]
      blockers = [(delta * a - anchor, delta * b - anchor) for (a, b) in blocker_template]

      if schedule == 'after':
        S = make_copies(T, offs) + blockers
      elif schedule == 'before':
        S = list(blockers) + make_copies(T, offs)
      elif schedule == 'interleaved':
        # add one copy then one blocker, etc., coupling copies and blockers tightly
        S = []
        for idx, start in enumerate(offs):
          S += make_copies(T, (start,))
          if idx < len(blockers):
            S.append(blockers[idx])
      else:  # 'split': half the copies, then blockers, then remaining copies
        h = max(1, len(offs) // 2)
        first = offs[:h]
        second = offs[h:]
        S = make_copies(T, first) + list(blockers) + make_copies(T, second)

      T = S
    return T

  # --- Enumerate candidates and select the best by ratio ---
  offsets_options = [
    (2, 6, 10, 14),         # baseline
    (1, 5, 9, 13),          # shifted pattern
    (3, 7, 11, 15),         # alternate stagger
    (2.5, 6.5, 10.5, 14.5), # half-step stagger to break symmetry
  ]
  translations = ['left', 'center']
  schedules = ['after', 'before', 'split', 'interleaved']
  extras = [False, True]
  depths = [max(2, iterations - 1), iterations]  # try one less and the given

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
            # build candidate (may be somewhat large; skip if blow-up occurs)
            T = build_candidate(k, schedule=sch, extra_first=ex, translation=trans, offsets=offs)
            if not T:
              continue
            # skip huge instances to keep search practical
            if len(T) > 4000:
              continue
            # normalize before testing so candidates are comparable
            Tn = _normalize_grid(T)
            om = clique_number(Tn)
            if om == 0:
              continue
            cols = firstfit_colors(Tn)
            ratio = cols / om
            n = len(Tn)
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
              best_T = Tn
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

  # Return the best found normalized intervals
  return best_T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()