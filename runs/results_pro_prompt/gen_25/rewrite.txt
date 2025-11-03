# EVOLVE-BLOCK-START

def construct_intervals():
  """
  Construct a sequence of open intervals presented to FirstFit.
  Strategy: depth-5 fractal/tile backbone with multiple anchors and
  cap placement variants under a tight candidate set. We return the
  best candidate's own arrival order to maximize FirstFit/OPT ratio.

  Returns:
    intervals: list of tuples (l, r) representing open intervals (l, r)
  """
  from typing import List, Tuple

  Interval = Tuple[float, float]

  # -------------------- Core utilities --------------------

  def clique_number(intervals: List[Interval]) -> int:
    # Sweep-line; ends before starts on ties (open intervals)
    events = []
    for (l, r) in intervals:
      events.append((l, 1))
      events.append((r, -1))
    events.sort(key=lambda x: (x[0], x[1]))
    active = 0
    best = 0
    for _, delta in events:
      active += delta
      if active > best:
        best = active
    return best

  def firstfit_color_count(intervals: List[Interval]) -> int:
    # FirstFit using color bins (each bin keeps its intervals)
    colors: List[List[Interval]] = []
    for (l, r) in intervals:
      placed = False
      # try to place in the first non-conflicting color
      for col in colors:
        conflict = False
        # open intervals overlap iff cl < r and l < cr
        for (cl, cr) in col:
          if cl < r and l < cr:
            conflict = True
            break
        if not conflict:
          col.append((l, r))
          placed = True
          break
      if not placed:
        colors.append([(l, r)])
    return len(colors)

  # -------------------- Fractal/Tiles builder --------------------

  def build_fractal(depth: int,
                    starts: Tuple[int, ...],
                    caps: List[Tuple[int, int]],
                    interleave: bool = False,
                    caps_before: bool = False) -> List[Interval]:
    """
    Recursive fractal-like builder:
    - starts: tuple of integer offsets for block copies per level
    - caps: list of (a,b) base intervals appended each level (scaled by delta)
    - interleave: copy order variant (kept False here to preserve canonical order)
    - caps_before: present caps before copies at each level or after
    """
    T: List[Interval] = [(0.0, 1.0)]
    for _ in range(depth):
      lo = min(l for l, r in T)
      hi = max(r for l, r in T)
      delta = hi - lo
      S: List[Interval] = []
      if caps_before:
        for a, b in caps:
          S.append((delta * a, delta * b))
      if interleave:
        for i in range(len(T)):
          li, ri = T[i]
          for s in starts:
            S.append((delta * s + li - lo, delta * s + ri - lo))
      else:
        for s in starts:
          for (l, r) in T:
            S.append((delta * s + l - lo, delta * s + r - lo))
      if not caps_before:
        for a, b in caps:
          S.append((delta * a, delta * b))
      T = S
    # Normalize so leftmost is 0.0
    lo = min(l for l, r in T)
    return [(l - lo, r - lo) for l, r in T]

  # -------------------- Candidate construction --------------------

  # Baseline caps from literature (Fig. 4 of the cited paper)
  baseline_caps = [(1, 5), (12, 16), (4, 9), (8, 13)]

  # Build a compact set of strong candidates at depth 5
  depth = 5  # gives 2388 intervals for 4 starts and 4 caps; within a reasonable bound for evaluator
  starts_list = [
    (2, 6, 10, 14),
    (3, 7, 11, 15),
    (4, 8, 12, 16),
  ]
  caps_before_opts = [False, True]

  candidates: List[List[Interval]] = []
  for starts in starts_list:
    for caps_before in caps_before_opts:
      seq = build_fractal(depth, starts, baseline_caps, interleave=False, caps_before=caps_before)
      candidates.append(seq)

  # Evaluate candidates and pick the one maximizing FirstFit/OPT
  best_seq = None
  best_ratio = -1.0

  for seq in candidates:
    omega = clique_number(seq)
    if omega <= 0:
      continue
    colors = firstfit_color_count(seq)
    ratio = colors / omega
    if ratio > best_ratio + 1e-12 or (abs(ratio - best_ratio) < 1e-12 and (best_seq is None or len(seq) < len(best_seq))):
      best_ratio = ratio
      best_seq = seq

  # Fallback (should not happen): depth-4 baseline
  if best_seq is None:
    best_seq = build_fractal(4, (2, 6, 10, 14), baseline_caps, interleave=False, caps_before=False)

  return best_seq

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()