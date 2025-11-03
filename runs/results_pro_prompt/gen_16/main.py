# EVOLVE-BLOCK-START

def construct_intervals():
  """
  Construct a sequence of intervals of real line,
  in the order in which they are presented to FirstFit,
  so that it maximizes the number of colors used by FirstFit
  divided by the maximum number of intervals that cover a single point

  The initial implementation uses the construction from
  Figure 4 in https://arxiv.org/abs/1506.00192

  Returns:
    intervals: list of tuples, each tuple (l, r) represents an open interval from l to r
  """

  import random
  random.seed(0)

  def firstfit_color_count(intervals):
    """Simulate FirstFit on the given arrival-ordered intervals."""
    colors = []
    for (l, r) in intervals:
      placed = False
      for col in colors:
        conflict = False
        for (cl, cr) in col:
          # open intervals overlap iff cl < r and l < cr
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

  def clique_number(intervals):
    """Compute the clique number (ω) by sweep-line. For open intervals we
    process ends before starts at the same coordinate so [a,b) and [b,c)
    don't overlap."""
    events = []
    for (l, r) in intervals:
      events.append((l, 1))   # start
      events.append((r, -1))  # end
    events.sort(key=lambda x: (x[0], x[1]))  # end (-1) comes before start (+1) at same coord
    active = 0
    best = 0
    for _, delta in events:
      active += delta
      if active > best:
        best = active
    return best

  def generate_recursive(starts, caps, depth, interleave=False, caps_before=False):
    """Generate a recursive construction parameterized by:
       - starts: sequence of integer offsets for block copies per level
       - caps: list of (a,b) base intervals appended each level (scaled by delta)
       - depth: recursion depth
       - interleave: if True, use a round-robin copy order (increases cross-copy interference)
       - caps_before: if True, present caps first at each level to pre-block low colors
    """
    T = [(0.0, 1.0)]
    for _ in range(depth):
      lo = min(l for l, r in T)
      hi = max(r for l, r in T)
      delta = hi - lo
      S = []
      if caps_before:
        # Place caps first this round to occupy low colors early
        for a, b in caps:
          S.append((delta * a, delta * b))
      if interleave:
        # Round-robin copy to increase color interactions across blocks
        for i in range(len(T)):
          li, ri = T[i]
          for start in starts:
            S.append((delta * start + li - lo, delta * start + ri - lo))
      else:
        for start in starts:
          S += [(delta * start + l - lo, delta * start + r - lo) for l, r in T]
      if not caps_before:
        # Append caps after block copies (original behavior)
        for a, b in caps:
          S.append((delta * a, delta * b))
      T = S
    return T

  def hybrid_blockwise(seq):
    """Block-wise hybrid ordering to exploit local fractal structure."""
    n = len(seq)
    if n < 8:
      return list(seq)
    k = 4
    b = n // k
    parts = [list(seq[i * b:(i + 1) * b]) for i in range(k - 1)]
    parts.append(list(seq[(k - 1) * b:]))
    out = []
    # block 0: keep original local structure
    out += parts[0]
    # block 1: left-first
    if len(parts) > 1:
      out += sorted(parts[1], key=lambda x: (x[0], x[1]))
    # block 2: reversed
    if len(parts) > 2:
      out += list(reversed(parts[2]))
    # block 3: short-first
    if len(parts) > 3:
      out += sorted(parts[3], key=lambda x: (x[1] - x[0], x[0]))
    return out

  # baseline original construction (kept as fallback)
  def original():
    T = [(0, 1)]
    for _ in range(4):
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
    return T

  # Search for an improved variant (time/size-bounded)
  best_seq = None
  best_score = -1.0
  best_meta = None
  max_intervals = 2200

  starts_options = [
    (2, 6, 10, 14),
    (3, 7, 11, 15),
    (4, 8, 12, 16),
  ]
  depth_options = [3, 4]
  trials_per_combo = 6

  # Start with the baseline as a safe candidate (evaluate multiple orders)
  baseline = original()
  baseline_omega = clique_number(baseline)
  if baseline_omega > 0:
    orders = []
    base = list(baseline)
    orders.append(("fractal", base))
    orders.append(("reversed", list(reversed(base))))
    orders.append(("left_first", sorted(base, key=lambda x: (x[0], x[1]))))
    orders.append(("right_first", sorted(base, key=lambda x: (-x[0], -x[1]))))
    orders.append(("short_first", sorted(base, key=lambda x: (x[1] - x[0], x[0]))))
    orders.append(("long_first", sorted(base, key=lambda x: (-(x[1] - x[0]), x[0]))))
    orders.append(("hybrid_blockwise", hybrid_blockwise(base)))
    for order_name, seq0 in orders:
      colors0 = firstfit_color_count(seq0)
      score0 = colors0 / baseline_omega
      if score0 > best_score:
        best_seq = seq0
        best_score = score0
        best_meta = ("baseline order=%s" % order_name, colors0, baseline_omega, len(seq0))

  # Randomized exploration of caps and interleaving (deterministic seed)
  for starts in starts_options:
    for depth in depth_options:
      # quick estimate of explosion: skip patterns that will definitely exceed budget
      p_est = max(4, len(starts))
      est = 1
      for _ in range(depth):
        est = len(starts) * est + p_est
        if est > max_intervals:
          break
      if est > max_intervals:
        continue

      min_s = min(starts)
      max_s = max(starts)

      # heuristic cap set derived from the starts sequence
      cap_heur = []
      for i in range(len(starts)):
        if i == 0:
          a = starts[0] - 1
          b = starts[0] + 3
        elif i == len(starts) - 1:
          a = starts[i - 1] + 2
          b = starts[i] + 2
        else:
          a = starts[i - 1] + 2
          b = starts[i] + 3
        cap_heur.append((a, b))

      candidates = [cap_heur]
      # add random cap-sets to explore the local neighborhood
      for _ in range(trials_per_combo):
        caps = []
        for _i in range(max(4, len(starts))):
          a = random.randint(min_s - 2, max_s + 1)
          b = a + random.randint(3, 6)
          caps.append((a, b))
        candidates.append(caps)

      for interleave in (False, True):
        for caps_before in (False, True):
          for caps in candidates:
            seq = generate_recursive(starts, caps, depth, interleave=interleave, caps_before=caps_before)
            if len(seq) == 0 or len(seq) > max_intervals:
              continue
            omega = clique_number(seq)
            if omega <= 0:
              continue

            # Evaluate several principled arrival orders
            base = list(seq)
            orderings = [
              ("fractal", base),
              ("reversed", list(reversed(base))),
              ("left_first", sorted(base, key=lambda x: (x[0], x[1]))),
              ("right_first", sorted(base, key=lambda x: (-x[0], -x[1]))),
              ("short_first", sorted(base, key=lambda x: (x[1] - x[0], x[0]))),
              ("long_first", sorted(base, key=lambda x: (-(x[1] - x[0]), x[0]))),
              ("hybrid_blockwise", hybrid_blockwise(base)),
            ]

            for order_name, candidate_seq in orderings:
              colors = firstfit_color_count(candidate_seq)
              score = colors / omega
              if score > best_score or (abs(score - best_score) < 1e-9 and len(candidate_seq) < len(best_seq)):
                best_score = score
                best_seq = candidate_seq
                best_meta = ("starts=%s depth=%d interleave=%s caps_before=%s caps=%s order=%s" % (starts, depth, interleave, caps_before, caps, order_name), colors, omega, len(candidate_seq))

  if best_seq is None:
    best_seq = original()
  return best_seq

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()