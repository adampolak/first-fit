# EVOLVE-BLOCK-START

def construct_intervals(iterations=4, extra_first=True):
  """
  Construct a sequence of intervals of real line,
  in the order in which they are presented to FirstFit,
  so that it maximizes the number of colors used by FirstFit
  divided by the maximum number of intervals that cover a single point

  This variant:
    - rotates copy offsets across four staggerings (A/B/C/D) over levels;
    - rotates among three blocker templates per level;
    - alternates anchoring: even levels left-anchored, odd levels center-anchored;
    - optionally adds an extra fifth copy on the first iteration to
      pollute early small colors (extra_first=True);
    - normalizes endpoints to a compact integer grid before returning.

  Arguments:
    iterations: number of recursive expansion steps (default 4)
    extra_first: whether to add a fifth copy at the first iteration (default True)

  Returns:
    intervals: list of tuples, each tuple (l, r) represents an open interval from l to r
  """

  # Four offset patterns (A/B/C/D) to cycle deterministically across levels
  offsets_patterns = [
    (2, 6, 10, 14),   # A: canonical
    (1, 5, 9, 13),    # B: shifted
    (3, 7, 11, 15),   # C: shifted
    (0, 4, 8, 12),    # D: shifted
  ]
  # Three blocker templates to rotate per level
  blocker_patterns = [
    ((1, 5), (12, 16), (4, 9), (8, 13)),   # canonical
    ((0, 4), (6, 10), (8, 12), (14, 18)),  # variant B
    ((2, 6), (4, 8), (10, 14), (12, 16)),  # variant C
  ]

  T = [(0.0, 1.0)]
  # Each iteration replaces T by multiple translated/scaled copies plus connectors.
  for lvl in range(iterations):
    lo = min(l for l, r in T)
    hi = max(r for l, r in T)
    delta = hi - lo
    center = (lo + hi) / 2.0
    S = []

    # pick offsets and blockers for this level (rotate patterns)
    base_offsets = offsets_patterns[lvl % len(offsets_patterns)]
    blockers = blocker_patterns[lvl % len(blocker_patterns)]
    if extra_first and lvl == 0:
      # extend first-level offsets with an extra copy to increase early coupling
      offs = tuple(list(base_offsets) + [18])
    else:
      offs = base_offsets

    # Choose anchoring rule: even levels left-anchored, odd levels center-anchored
    center_anchor = (lvl % 2 == 1)

    # Place translated copies
    for start in offs:
      if center_anchor:
        offset = delta * start - center
      else:
        offset = delta * start - lo
      S += [(l + offset, r + offset) for (l, r) in T]

    # Add connecting long intervals (blockers), anchored consistently with this level
    if center_anchor:
      S += [(delta * a - center, delta * b - center) for (a, b) in blockers]
    else:
      S += [(delta * a, delta * b) for (a, b) in blockers]

    T = S

  # Normalize endpoints to a compact integer grid while preserving order.
  # Map each unique endpoint to an increasing even integer (to ensure
  # positive lengths and avoid degeneracy).
  endpoints = sorted(set([x for seg in T for x in seg]))
  coord = {}
  cur = 0
  for e in endpoints:
    coord[e] = cur
    cur += 2
  normalized = [(coord[l], coord[r]) for (l, r) in T]

  return normalized

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