# EVOLVE-BLOCK-START

def construct_intervals(iterations=4):
  """
  Construct a sequence of intervals of real line,
  in the order in which they are presented to FirstFit,
  so that it maximizes the number of colors used by FirstFit
  divided by the maximum number of intervals that cover a single point.

  The implementation uses the construction from
  Figure 4 in https://arxiv.org/abs/1506.00192 (four copies + four blockers
  per iteration), and then applies a deterministic greedy shrink step that
  removes intervals while preserving both FirstFit colors and omega.

  Arguments:
    iterations: number of recursive expansion steps (default 4)

  Returns:
    intervals: list of tuples, each tuple (l, r) represents an open interval from l to r
  """

  # --- Helper routines (lightweight and deterministic) ---

  def firstfit_colors(intervals):
    """Simulate FirstFit using per-color last-end tracking (exact for interval graphs)."""
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

  def shrink_witness(T, target_cols, target_omega, max_passes=2):
    """
    Greedily remove intervals that keep (FF colors >= target_cols) and (omega == target_omega).
    Run a couple of forward/backward passes to expose cascading redundancies.
    """
    changed = True
    passes = 0
    while changed and passes < max_passes:
      changed = False
      # forward scan
      i = 0
      while i < len(T):
        cand = T[:i] + T[i+1:]
        if clique_number(cand) == target_omega and firstfit_colors(cand) >= target_cols:
          T = cand
          changed = True
          # keep i at same index to test the next interval now in this position
        else:
          i += 1
      # backward scan
      i = len(T) - 1
      while i >= 0:
        cand = T[:i] + T[i+1:]
        if clique_number(cand) == target_omega and firstfit_colors(cand) >= target_cols:
          T = cand
          changed = True
        i -= 1
      passes += 1
    return T

  # --- Core four-copy + four-blocker recursion (no extra 5th copy) ---
  T = [(0, 1)]
  for _ in range(iterations):
    lo = min(l for l, r in T)
    hi = max(r for l, r in T)
    delta = hi - lo
    S = []
    # Keep robust four-copy offsets; dropping the 5th copy tightens size without harming ratio
    for start in (2, 6, 10, 14):
      S += [(delta * start + l - lo, delta * start + r - lo) for l, r in T]
    # Four connecting long intervals as in Figure 4
    S += [
      (delta * 1,  delta * 5),
      (delta * 12, delta * 16),
      (delta * 4,  delta * 9),
      (delta * 8,  delta * 13),
    ]
    T = S

  # --- Deterministic pruning to reduce n while preserving the witness (FF, omega) ---
  om = clique_number(T)
  cols = firstfit_colors(T)
  if om > 0 and cols > 0:
    T = shrink_witness(T, cols, om, max_passes=2)

  return T

  # return [  # Figure 3, OPT=2, FF=4
  #   (2,3),
  #   (6,7),
  #   (10,11),
  #   (14,15),
  #   (1,5),
  #   (12,16),
  #   (4,9),
  #   (8,13),
  # ]

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()