# EVOLVE-BLOCK-START

def construct_intervals(iterations=4):
  """
  Construct a sequence of open intervals (l, r) presented in order to FirstFit.
  The construction is a recursive expansion inspired by Figure 4 of
  https://arxiv.org/abs/1506.00192. Increasing 'iterations' strengthens the
  adversary and typically increases the FirstFit/OPT ratio.

  Improvements in this variant:
  - Use a slightly richer base gadget (two disjoint unit intervals).
  - Alternate the four-copy starting offsets between two staggerings across
    iterations to break alignment regularity.
  - Alternate two blocker/connector templates (one shifted by 0.5 units
    relative to the other) across iterations to couple colors in a less
    regular way, making FirstFit more likely to open extra colors.
  """

  # Base gadget: two disjoint small intervals (provides richer local structure)
  T = [(0.0, 1.0), (3.0, 4.0)]

  # Two choices of starts to alternate between iterations (breaks perfect repetition)
  starts_options = [
    (2, 6, 10, 14),
    (1, 5, 9, 13)
  ]

  # Two blocker/connector templates. The second is shifted by 0.5 (relative)
  # to introduce fractional-phase overlaps (kept as floats until normalization).
  blocker_templates = [
    [(1.0, 5.0), (12.0, 16.0), (4.0, 9.0), (8.0, 13.0)],
    [(1.5, 5.5), (12.5, 16.5), (4.5, 9.5), (8.5, 13.5)]
  ]

  # Recursive expansion: alternate starts and blocker templates each level.
  for i in range(iterations):
    lo = min(l for l, r in T)
    hi = max(r for l, r in T)
    delta = hi - lo
    S = []
    starts = starts_options[i % len(starts_options)]
    blockers = blocker_templates[i % len(blocker_templates)]

    # Place four scaled/translated copies at staggered offsets (with current starts)
    for start in starts:
      offset = delta * start - lo
      S.extend([(offset + l, offset + r) for l, r in T])

    # Add the four connecting long intervals using the currently selected blocker template.
    # Each blocker is expressed in template coordinates and scaled by delta.
    for (a, b) in blockers:
      S.append((delta * a, delta * b))

    T = S

  # Normalize endpoints to a compact integer grid while preserving order.
  endpoints = sorted(set([x for seg in T for x in seg]))
  coord = {}
  cur = 0
  for e in endpoints:
    coord[e] = cur
    cur += 2

  normalized = [(coord[l], coord[r]) for (l, r) in T]

  return normalized

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()