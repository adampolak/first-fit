# EVOLVE-BLOCK-START

def construct_intervals():
  """
  Construct a sequence of intervals of real line,
  in the order in which they are presented to FirstFit,
  so that it maximizes the number of colors used by FirstFit
  divided by the maximum number of intervals that cover a single point.

  Strategy:
  - Diversified fractal/tile recursive construction with multiple anchor families.
  - Budget-aware pruning on recursion depth and cap multiplicity.
  - Multi-order screening (deterministic + random shuffles).
  - Bitset-based adversarial ordering refinement to force new FirstFit colors.
  - Internal cross-verification for ω with two independent methods.

  Returns:
    intervals: list of tuples, each tuple (l, r) represents an open interval [l, r)
  """
  import random
  from bisect import bisect_left
  from typing import List, Tuple

  random.seed(0)

  Interval = Tuple[float, float]

  # -------------------- FirstFit (optimized) --------------------

  def firstfit_color_count(intervals: List[Interval]) -> int:
    """
    Simulate FirstFit on the given arrival-ordered intervals.
    Optimization: maintain each color as a list sorted by left endpoint,
    and check only immediate neighbors for conflicts (open-interval semantics).
    """
    colors: List[List[Interval]] = []

    for (l, r) in intervals:
      placed = False
      for col in colors:
        # Find insertion position by left endpoint
        pos = bisect_left(col, (l, r))
        # Check previous interval in color
        if pos > 0 and col[pos - 1][1] > l:  # open: overlap if prev.right > l
          continue
        # Check next interval in color
        if pos < len(col) and col[pos][0] < r:  # open: overlap if next.left < r
          continue
        # Can place here
        col.insert(pos, (l, r))
        placed = True
        break
      if not placed:
        colors.append([(l, r)])
    return len(colors)

  # -------------------- Clique number ω --------------------

  def clique_number_sweepline(intervals: List[Interval]) -> int:
    """
    Sweep-line for open intervals (process ends before starts on ties).
    """
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

  def clique_number_sampling(intervals: List[Interval]) -> int:
    """
    Independent ω check: sample midpoints between consecutive unique endpoints.
    Max overlap occurs on some open segment between distinct endpoints.
    """
    pts = set()
    for (l, r) in intervals:
      pts.add(l); pts.add(r)
    xs = sorted(pts)
    if len(xs) <= 1:
      return len(intervals)
    # Sample midpoints between consecutive endpoints
    samples = [0.5 * (xs[i] + xs[i + 1]) for i in range(len(xs) - 1)]
    best = 0
    # Two-pointer approach: pre-sort by l and r
    by_l = sorted(intervals, key=lambda x: x[0])
    by_r = sorted(intervals, key=lambda x: x[1])
    L = R = 0
    active = 0
    # Sweep over samples; maintain active intervals respecting open semantics
    # We'll advance pointers to include starts < s and remove ends <= s
    for s in samples:
      while L < len(by_l) and by_l[L][0] < s:
        # Interval starts before s; it may still be active if its end > s
        active += 1
        L += 1
      while R < len(by_r) and by_r[R][1] <= s:
        # Interval ended at or before s (open intervals exclude endpoint)
        active -= 1
        R += 1
      if active > best:
        best = active
    return best

  def clique_number(intervals: List[Interval]) -> int:
    """
    Cross-verify ω. Prefer sweep-line; if mismatch occurs, take the maximum.
    """
    a = clique_number_sweepline(intervals)
    b = clique_number_sampling(intervals)
    return max(a, b)

  # -------------------- Orders --------------------

  def order_identity(seq: List[Interval]) -> List[Interval]:
    return list(seq)

  def order_reversed(seq: List[Interval]) -> List[Interval]:
    return list(reversed(seq))

  def order_left_first(seq: List[Interval]) -> List[Interval]:
    return sorted(seq, key=lambda x: (x[0], x[1]))

  def order_right_first(seq: List[Interval]) -> List[Interval]:
    return sorted(seq, key=lambda x: (-x[0], -x[1]))

  def order_short_first(seq: List[Interval]) -> List[Interval]:
    return sorted(seq, key=lambda x: ((x[1] - x[0]), x[0]))

  def order_long_first(seq: List[Interval]) -> List[Interval]:
    return sorted(seq, key=lambda x: (-(x[1] - x[0]), x[0]))

  def order_by_right_endpoint(seq: List[Interval]) -> List[Interval]:
    return sorted(seq, key=lambda x: (x[1], x[0]))

  def order_hybrid_blockwise(seq: List[Interval]) -> List[Interval]:
    n = len(seq)
    if n < 16:
      return list(seq)
    k = 4
    b = n // k
    parts = [list(seq[i * b:(i + 1) * b]) for i in range(k - 1)]
    parts.append(list(seq[(k - 1) * b:]))
    out = []
    out += parts[0]  # original local
    if len(parts) > 1:
      out += order_left_first(parts[1])
    if len(parts) > 2:
      out += order_reversed(parts[2])
    if len(parts) > 3:
      out += order_short_first(parts[3])
    return out

  def order_alt_short_long(seq: List[Interval]) -> List[Interval]:
    n = len(seq)
    if n <= 2:
      return list(seq)
    by_len = sorted(seq, key=lambda x: ((x[1] - x[0]), x[0]))
    i, j = 0, n - 1
    out = []
    take_small = True
    while i <= j:
      if take_small:
        out.append(by_len[i]); i += 1
      else:
        out.append(by_len[j]); j -= 1
      take_small = not take_small
    return out

  ORDERERS = [
    ("identity", order_identity),
    ("reversed", order_reversed),
    ("left_first", order_left_first),
    ("right_first", order_right_first),
    ("short_first", order_short_first),
    ("long_first", order_long_first),
    ("right_endpoint", order_by_right_endpoint),
    ("hybrid_blockwise", order_hybrid_blockwise),
    ("alt_short_long", order_alt_short_long),
  ]

  # -------------------- Adversarial ordering --------------------

  def adversarial_order(seq: List[Interval], seed: int = 17) -> List[Interval]:
    """
    Adaptive adversary that greedily picks the next interval to force a new color
    in FirstFit if possible. Uses neighbor bitsets to track conflicts against
    current color classes.
    """
    n = len(seq)
    if n == 0:
      return []

    # Precompute overlaps bitsets
    neighbors = [0] * n
    for i in range(n):
      li, ri = seq[i]
      for j in range(i + 1, n):
        lj, rj = seq[j]
        if li < rj and lj < ri:
          neighbors[i] |= (1 << j)
          neighbors[j] |= (1 << i)

    unpicked = set(range(n))
    color_sets = []  # list of bitsets for each color
    order_idx = []
    rng = random.Random(seed)

    def blocks_all_colors(i: int) -> bool:
      ni = neighbors[i]
      for cmask in color_sets:
        if (ni & cmask) == 0:
          return False
      return True

    def blocked_colors_count(i: int) -> int:
      ni = neighbors[i]
      cnt = 0
      for cmask in color_sets:
        if (ni & cmask) != 0:
          cnt += 1
      return cnt

    def future_coverage(i: int) -> int:
      # Count overlaps with still-unpicked intervals
      ni = neighbors[i]
      mask_unpicked = 0
      for j in unpicked:
        mask_unpicked |= (1 << j)
      if i in unpicked:
        mask_unpicked &= ~(1 << i)
      return (ni & mask_unpicked).bit_count()

    def assign_color(i: int):
      ni = neighbors[i]
      for c in range(len(color_sets)):
        if (ni & color_sets[c]) == 0:
          color_sets[c] |= (1 << i)
          return
      # need a new color
      color_sets.append(1 << i)

    while unpicked:
      used = len(color_sets)
      best = None
      best_key = None
      # small random sample if very large, else all
      picks = list(unpicked) if len(unpicked) < 200 else rng.sample(list(unpicked), 200)
      for i in picks:
        blk_all = 1 if blocks_all_colors(i) else 0
        blk_cnt = blocked_colors_count(i)
        cov = future_coverage(i)
        length = seq[i][1] - seq[i][0]
        key = (blk_all, blk_cnt, cov, -length)
        if best_key is None or key > best_key:
          best_key = key
          best = i
      order_idx.append(best)
      assign_color(best)
      unpicked.remove(best)

    return [seq[i] for i in order_idx]

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
    - interleave: copy order variant
    - caps_before: present caps before copies at each level
    """
    T = [(0.0, 1.0)]
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
    return T  # keep absolute scale (no normalization) for diversity

  def cap_templates_for_starts(starts: Tuple[int, ...]) -> List[List[Tuple[int, int]]]:
    """
    Produce several cap sets for the given starts. Aim to bridge adjacent
    and non-adjacent blocks to increase FirstFit colors without exploding ω.
    """
    s = list(starts)
    t: List[List[Tuple[int, int]]] = []
    m = len(s)
    if m >= 4:
      # Canonical set (paper-inspired)
      t.append([(s[0] - 1, s[0] + 3),
                (s[1] + 2, s[2] + 2),
                (s[2] - 1, s[3] + 1),
                (s[1] - 2, s[2] + 3)])
      # Balanced bridges
      t.append([(s[0], s[1] + 2),
                (s[1] + 1, s[2] + 3),
                (s[2], s[3] + 2),
                (s[0] + 2, s[2] + 1)])
      # Wider caps
      t.append([(s[0] - 1, s[1] + 3),
                (s[1] - 1, s[2] + 3),
                (s[2] - 1, s[3] + 3),
                (s[0] + 1, s[3] + 1)])
      # Sparse caps
      t.append([(s[0] - 1, s[1] + 1),
                (s[2] - 1, s[3] + 1),
                (s[1], s[2] + 2),
                (s[0] + 2, s[3])])
    else:
      # Simple adjacent bridges replicated
      for i in range(m - 1):
        t.append([(s[i], s[i + 1] + 2)] * max(4, m))
    return t

  # -------------------- Search pipeline --------------------

  def evaluate_orders(seq: List[Interval], orderers, rnd_orders=4, rnd_seed=2025):
    """
    Evaluate multiple deterministic orders plus a few seeded random shuffles.
    Return:
      best_colors, best_ordered_seq, mean_ratio_over_orders (w.r.t. fixed ω)
    """
    omega = clique_number(seq)
    if omega <= 0:
      return (0, list(seq), 0.0, omega)

    best_colors = -1
    best_seq = None

    # Deterministic orders
    colors_list = []
    for name, fn in orderers:
      ord_seq = fn(seq)
      c = firstfit_color_count(ord_seq)
      colors_list.append(c)
      if c > best_colors:
        best_colors = c
        best_seq = ord_seq

    # Seeded random shuffles
    rng = random.Random(rnd_seed + len(seq))
    for k in range(rnd_orders):
      tmp = list(seq)
      rng.shuffle(tmp)
      c = firstfit_color_count(tmp)
      colors_list.append(c)
      if c > best_colors:
        best_colors = c
        best_seq = tmp

    mean_ratio = sum(c / omega for c in colors_list) / len(colors_list)
    return (best_colors, best_seq, mean_ratio, omega)

  # Limits and configuration
  max_intervals = 740  # strict budget
  depths = [3, 4, 5]
  starts_families = [
    (2, 6, 10, 14),
    (3, 7, 11, 15),
    (4, 8, 12, 16),
  ]
  bool_flags = [(False, False), (True, False), (False, True), (True, True)]

  # Baseline caps (literature)
  baseline_caps = [(1, 5), (12, 16), (4, 9), (8, 13)]

  # Candidate generation
  candidates: List[List[Interval]] = []

  # Always include three anchored baselines at depth 4
  for st in starts_families:
    seq = build_fractal(4, st, baseline_caps, interleave=False, caps_before=False)
    if len(seq) <= max_intervals:
      candidates.append(seq)

  # Controlled exploration with templates + jitter
  rng = random.Random(1337)
  for starts in starts_families:
    templates = cap_templates_for_starts(starts)
    jittered = []
    for tpl in templates:
      for _ in range(2):  # two jitter variants per template
        caps = []
        for (a, b) in tpl:
          da = rng.choice([-2, -1, 0, 1, 2])
          db = rng.choice([-1, 0, 1, 2, 3])
          aa, bb = a + da, b + db
          if bb - aa < 3:
            bb = aa + 3
          caps.append((aa, bb))
        jittered.append(caps)
    caps_sets = templates + jittered

    for depth in depths:
      for caps in caps_sets:
        for interleave, caps_before in bool_flags:
          # Budget-aware size estimate
          est = 1
          feasible = True
          for _ in range(depth):
            est = len(starts) * est + len(caps)
            if est > max_intervals:
              feasible = False
              break
          if not feasible or est == 0:
            continue
          seq = build_fractal(depth, starts, caps, interleave=interleave, caps_before=caps_before)
          if 0 < len(seq) <= max_intervals:
            candidates.append(seq)

  # Remove duplicates by string fingerprint (lightweight)
  seen = set()
  unique_candidates = []
  for seq in candidates:
    key = (len(seq), round(min(l for l, r in seq), 6), round(max(r for l, r in seq), 6))
    if key not in seen:
      seen.add(key)
      unique_candidates.append(seq)

  # Screen candidates with multi-order evaluation; keep top by a mix of best ratio and robustness
  scored = []
  for seq in unique_candidates:
    best_colors, best_seq_ord, mean_ratio, omega = evaluate_orders(seq, ORDERERS, rnd_orders=4)
    if omega <= 0:
      continue
    best_ratio = best_colors / omega
    # Score: prioritize best_ratio; tiebreak by mean_ratio, then shorter length
    scored.append((best_ratio, mean_ratio, len(seq), seq, best_seq_ord, omega, best_colors))

  if not scored:
    # Fallback to depth-4 baseline
    return build_fractal(4, (2, 6, 10, 14), baseline_caps, interleave=False, caps_before=False)

  scored.sort(key=lambda x: (-x[0], -x[1], x[2]))
  top_k = min(6, len(scored))
  top = scored[:top_k]

  # Adversarial refinement on top candidates
  best_overall_ratio = -1.0
  best_overall_seq = None
  for best_ratio, mean_ratio, nlen, seq_set, seq_ord, omega, colors_det in top:
    # Run adversarial order on raw set and on the best deterministic order
    seq_adv = adversarial_order(seq_set, seed=101)
    colors_adv = firstfit_color_count(seq_adv)
    ratio_adv = colors_adv / omega

    seq_adv2 = adversarial_order(seq_ord, seed=202)
    colors_adv2 = firstfit_color_count(seq_adv2)
    ratio_adv2 = colors_adv2 / omega

    # Also try combining adversarial with blockwise hybrid to reorder within blocks
    seq_hyb_adv = order_hybrid_blockwise(seq_adv)
    colors_hyb_adv = firstfit_color_count(seq_hyb_adv)
    ratio_hyb_adv = colors_hyb_adv / omega

    # Select best among deterministic and adversarial refinements
    candidates_final = [
      (best_ratio, seq_ord),
      (ratio_adv, seq_adv),
      (ratio_adv2, seq_adv2),
      (ratio_hyb_adv, seq_hyb_adv),
    ]
    for rr, seq_fin in candidates_final:
      if rr > best_overall_ratio or (abs(rr - best_overall_ratio) < 1e-12 and len(seq_fin) < len(best_overall_seq or [])):
        best_overall_ratio = rr
        best_overall_seq = seq_fin

  # Final safety: ensure ω is positive; if not, fallback
  if not best_overall_seq:
    best_overall_seq = build_fractal(4, (2, 6, 10, 14), baseline_caps, interleave=False, caps_before=False)

  # Internal verification (silent)
  _ff = firstfit_color_count(best_overall_seq)
  _opt = clique_number(best_overall_seq)
  # Return final sequence (arrival order)
  return best_overall_seq

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()