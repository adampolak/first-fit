# EVOLVE-BLOCK-START

def construct_intervals(iterations=4):
  """
  Construct a sequence of intervals of real line to maximize FirstFit/omega.
  Improvements over the previous implementation:
   - try both left-anchored and center-anchored translations for copies,
     and a richer set of alternative offset patterns to diversify geometry;
   - multiple blocker templates (A/B/C/D) to alter coupling between copies;
   - optional post-iteration cap intervals to couple colors across layers;
   - try deeper recursion depths (k up to iterations+2);
   - normalize endpoints to a compact even integer grid before returning.
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

  def build_candidate(k, schedule='after', extra_first=False, translation='left',
                      offsets=(2,6,10,14), blocker_template='A', caps='none'):
    """
    Build k-level recursive pattern using copies and blockers.
    translation in {'left','center'}:
      - 'left' anchors copies by left endpoint as before
      - 'center' aligns copies relative to the center of current T,
        producing a staggered geometry that can be adversarial to FirstFit
    offsets is a tuple of integer multipliers for the level translation.
    blocker_template in {'A','B','C','D'} selects the blocker geometry.
    caps in {'none','cap1','cap2'} injects optional post-iteration long intervals.
    """
    T = [(0.0, 1.0)]
    for i in range(k):
      lo = min(l for l, r in T)
      hi = max(r for l, r in T)
      delta = hi - lo
      center = (lo + hi) / 2.0

      # optionally include a fifth copy on the first iteration
      offs = tuple(list(offsets) + ([18] if extra_first and i == 0 else []))

      def make_copies(from_T, offsets_local):
        S = []
        for start in offsets_local:
          if translation == 'left':
            offset = delta * start - lo
          else:  # center translation: align by center
            offset = delta * start - center
          for (l, r) in from_T:
            S.append((l + offset, r + offset))
        return S

      # Blocker templates (scaled by delta)
      if blocker_template == 'A':
        blockers = [
          (delta * 1,  delta * 5),
          (delta * 12, delta * 16),
          (delta * 4,  delta * 9),
          (delta * 8,  delta * 13),
        ]
      elif blocker_template == 'B':
        # Slightly left-shifted to change overlap tensions
        blockers = [
          (delta * 0,  delta * 4),
          (delta * 11, delta * 15),
          (delta * 3,  delta * 8),
          (delta * 7,  delta * 12),
        ]
      elif blocker_template == 'C':
        # Slightly right-shifted variant
        blockers = [
          (delta * 2,  delta * 6),
          (delta * 13, delta * 17),
          (delta * 5,  delta * 10),
          (delta * 9,  delta * 14),
        ]
      else:  # 'D' aggressive cross coupling
        blockers = [
          (delta * 0,  delta * 9),
          (delta * 8,  delta * 17),
          (delta * 1,  delta * 12),
          (delta * 4,  delta * 13),
        ]

      if schedule == 'after':
        S = make_copies(T, offs) + blockers
      elif schedule == 'before':
        S = list(blockers) + make_copies(T, offs)
      else:  # 'split': half the copies, then blockers, then remaining copies
        h = max(1, len(offs) // 2)
        first = offs[:h]
        second = offs[h:]
        S = make_copies(T, first) + list(blockers) + make_copies(T, second)

      T = S

    # Optional caps injected after final iteration to couple colors across layers
    if caps != 'none' and T:
      lo = min(l for l, r in T)
      hi = max(r for l, r in T)
      d = hi - lo
      if caps == 'cap1':
        T = T + [
          (d * 2,  d * 7),
          (d * 6,  d * 11),
          (d * 10, d * 15),
        ]
      elif caps == 'cap2':
        T = T + [
          (d * 1,  d * 3),
          (d * 5,  d * 7),
          (d * 9,  d * 11),
          (d * 13, d * 15),
        ]
    return T

  # --- Enumerate a modest set of candidates and select the best by ratio ---
  depths = list({max(2, iterations - 1), iterations, iterations + 1, iterations + 2})
  schedules = ['after', 'before', 'split']
  extras = [False, True]
  translations = ['left', 'center']
  offsets_options = [
    (2, 6, 10, 14),  # baseline
    (1, 5, 9, 13),   # shifted pattern
    (3, 7, 11, 15),  # alternate stagger
    (0, 4, 8, 12),   # zero-anchored stagger
  ]
  blocker_templates = ['A', 'B', 'C', 'D']
  caps_options = ['none', 'cap1', 'cap2']

  best_T = None
  best_ratio = -1.0
  best_n = None
  best_cols = 0
  best_om = 0

  for k in depths:
    for sch in schedules:
      for ex in extras:
        for trans in translations:
          for offs in offsets_options:
            for bt in blocker_templates:
              for cap in caps_options:
                T = build_candidate(k, schedule=sch, extra_first=ex,
                                    translation=trans, offsets=offs,
                                    blocker_template=bt, caps=cap)
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
    return _normalize_grid(T)

  # Normalize the best candidate to a compact integer grid before returning.
  return _normalize_grid(best_T)

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()