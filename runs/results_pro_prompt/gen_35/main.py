# EVOLVE-BLOCK-START

def construct_intervals():
  """
  Construct a sequence of intervals of real line,
  in the order in which they are presented to FirstFit,
  so that it maximizes the number of colors used by FirstFit
  divided by the maximum number of intervals that cover a single point

  The implementation is a cross-over of two strategies:
  - fractal/tile recursive construction
  - caching-based evaluation with multiple anchor strategies
  Returns:
    intervals: list of tuples, each tuple (l, r) represents an open interval [l, r)
  """
  import random
  random.seed(0)
  from functools import lru_cache

  @lru_cache(maxsize=None)
  def eval_seq_cached(seq_tuple):
    """Return (colors_used, omega) for a given sequence tuple representation.

    Improvement: evaluate multiple canonical orderings and an adversarial
    ordering heuristic; return the maximum FirstFit color count found and the
    clique-number (omega). This steers the search toward sequences that are
    strong against adversarial arrival orders.
    """
    seq = list(seq_tuple)

    # Compute clique number (omega) via sweep-line for open intervals
    events = []
    for (l, r) in seq:
      events.append((l, 1))   # start
      events.append((r, -1))  # end
    events.sort(key=lambda x: (x[0], x[1]))
    active = 0
    best = 0
    for _, delta in events:
      active += delta
      if active > best:
        best = active
    omega = best

    # Local FirstFit implementation (order-sensitive)
    def local_firstfit(sequence):
      colors = []
      for (l, r) in sequence:
        placed = False
        for col in colors:
          conflict = False
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

    # Canonical orderers to probe
    def order_identity(s): return list(s)
    def order_reversed(s): return list(reversed(s))
    def order_left_first(s): return sorted(s, key=lambda x: (x[0], x[1]))
    def order_right_first(s): return sorted(s, key=lambda x: (-x[0], -x[1]))
    def order_short_first(s): return sorted(s, key=lambda x: (x[1] - x[0], x[0]))
    def order_long_first(s): return sorted(s, key=lambda x: (-(x[1] - x[0]), x[0]))
    def order_right_endpoint(s): return sorted(s, key=lambda x: (x[1], x[0]))

    def order_hybrid_blockwise(s):
      n = len(s)
      if n < 12:
        return list(s)
      k = 4
      bsize = n // k
      parts = [list(s[i * bsize:(i + 1) * bsize]) for i in range(k - 1)]
      parts.append(list(s[(k - 1) * bsize:]))
      out = []
      out += parts[0]
      if len(parts) > 1:
        out += order_left_first(parts[1])
      if len(parts) > 2:
        out += list(reversed(parts[2]))
      if len(parts) > 3:
        out += order_short_first(parts[3])
      return out

    ORDERERS = [
      ("identity", order_identity),
      ("reversed", order_reversed),
      ("left_first", order_left_first),
      ("right_first", order_right_first),
      ("short_first", order_short_first),
      ("long_first", order_long_first),
      ("right_endpoint", order_right_endpoint),
      ("hybrid_blockwise", order_hybrid_blockwise),
    ]

    # Evaluate canonical orders and keep the maximum FirstFit color count
    best_colors = -1
    best_order_name = None
    for name, ord_fn in ORDERERS:
      try:
        seq_ord = ord_fn(seq)
        colors = local_firstfit(seq_ord)
        if colors > best_colors:
          best_colors = colors
          best_order_name = name
      except Exception:
        # be robust to any unexpected ordering failure
        continue

    # Adversarial ordering heuristic (bitset-based). Tries to choose the next
    # arriving interval that forces a new color (if possible) or otherwise
    # maximizes blocked colors / future coverage.
    def adversarial_order_local(sequence):
      n = len(sequence)
      if n == 0:
        return []
      # Precompute neighbor bitsets (Python int used as bitset)
      neighbors = [0] * n
      for i in range(n):
        li, ri = sequence[i]
        mask = 0
        for j in range(i + 1, n):
          lj, rj = sequence[j]
          if li < rj and lj < ri:
            mask |= (1 << j)
            neighbors[j] |= (1 << i)
        neighbors[i] |= mask

      unplaced = set(range(n))
      color_sets = []  # list of bitmasks of indices per color
      order_indices = []

      while unplaced:
        best_i = None
        best_key = None
        # compose mask of unplaced
        mask_unplaced = 0
        for j in unplaced:
          mask_unplaced |= (1 << j)
        used_colors = len(color_sets)

        for i in list(unplaced):
          ni = neighbors[i]
          # count how many current colors are blocked by i
          block_cnt = 0
          for cmask in color_sets:
            if (ni & cmask) != 0:
              block_cnt += 1
          new_color = (block_cnt == used_colors)
          coverage = ((ni & mask_unplaced) & ~(1 << i)).bit_count()
          length = sequence[i][1] - sequence[i][0]
          # lexicographic preference: force new_color, more blocks, more coverage, then shorter
          key = (1 if new_color else 0, block_cnt, coverage, -length)
          if best_key is None or key > best_key:
            best_key = key
            best_i = i

        # assign best_i to a FirstFit color in our simulated state
        ni = neighbors[best_i]
        placed = False
        for c, cmask in enumerate(color_sets):
          if (ni & cmask) == 0:
            color_sets[c] |= (1 << best_i)
            placed = True
            break
        if not placed:
          color_sets.append(1 << best_i)
        order_indices.append(best_i)
        unplaced.remove(best_i)

      return [sequence[i] for i in order_indices]

    # Try adversarial ordering if sequence isn't huge (safeguard)
    try:
      seq_adv = adversarial_order_local(seq)
      colors_adv = local_firstfit(seq_adv)
      if colors_adv > best_colors:
        best_colors = colors_adv
        best_order_name = "adversarial"
    except Exception:
      # ignore adversarial failures; fall back to canonical best
      pass

    # If nothing found (shouldn't happen), fall back to trivial values
    if best_colors < 0:
      best_colors = local_firstfit(seq)

    return (best_colors, omega)

  def firstfit_color_count(intervals):
    """Simulate FirstFit on the given arrival-ordered intervals."""
    colors = []
    for (l, r) in intervals:
      placed = False
      for col in colors:
        conflict = False
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

  def clique_number(intervals):
    """Compute the clique number (ω) by sweep-line for open intervals."""
    events = []
    for (l, r) in intervals:
      events.append((l, 1))   # start
      events.append((r, -1))  # end
    events.sort(key=lambda x: (x[0], x[1]))
    active = 0
    best = 0
    for _, delta in events:
      active += delta
      if active > best:
        best = active
    return best

  def generate_recursive(starts, caps, depth, interleave=False, caps_before=False):
    """Recursive fractal-like construction with optional interleaving and cap placement ordering."""
    T = [(0.0, 1.0)]
    for _ in range(depth):
      lo = min(l for l, r in T)
      hi = max(r for l, r in T)
      delta = hi - lo
      S = []
      if caps_before:
        # place caps first to occupy some colors early
        for a, b in caps:
          S.append((delta * a, delta * b))
      if interleave:
        for i in range(len(T)):
          li, ri = T[i]
          for start in starts:
            S.append((delta * start + li - lo, delta * start + ri - lo))
      else:
        for start in starts:
          S += [(delta * start + l - lo, delta * start + r - lo) for l, r in T]
      if not caps_before:
        for a, b in caps:
          S.append((delta * a, delta * b))
      T = S
    return T

  # Baseline original construction (safe fallback)
  def original():
    T = [(0.0, 1.0)]
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

  best_seq = None
  best_score = -1.0
  best_meta = None
  max_intervals = 3200

  # Variation space
  starts_options = [
    (2, 6, 10, 14),
    (3, 7, 11, 15),
    (2, 6, 10, 14, 18),
    (2, 6, 10, 14, 18, 22),
    (3, 7, 11, 15, 19),
  ]
  depth_options = [3, 4, 5]
  trials_per_combo = 60

  # Baseline evaluation
  baseline = original()
  baseline_colors, baseline_omega = eval_seq_cached(tuple(baseline))
  if baseline_omega > 0:
    best_seq = baseline
    best_score = baseline_colors / baseline_omega
    best_meta = ("baseline", baseline_colors, baseline_omega, len(baseline))

  # Alt anchor variant (3,7,11,15)
  def original_alt():
    T = [(0.0, 1.0)]
    for _ in range(4):
      lo = min(l for l, r in T)
      hi = max(r for l, r in T)
      delta = hi - lo
      S = []
      for start in (3, 7, 11, 15):
        S += [(delta * start + l - lo, delta * start + r - lo) for l, r in T]
      S += [
        (delta * 1, delta * 5),
        (delta * 12, delta * 16),
        (delta * 4, delta * 9),
        (delta * 8, delta * 13)
      ]
      T = S
    return T

  alt = original_alt()
  alt_colors, alt_omega = eval_seq_cached(tuple(alt))
  if alt_omega > 0:
    alt_score = alt_colors / alt_omega
    if alt_score > best_score + 1e-12:
      best_seq = alt
      best_score = alt_score
      best_meta = ("alt_anchor=(3,7,11,15)", alt_colors, alt_omega, len(alt))

  # Exploration over combinations
  for starts in starts_options:
    for depth in depth_options:
      # rough estimate of produced size; prune if exceeds budget
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
            colors, omega = eval_seq_cached(tuple(seq))
            if omega <= 0:
              continue
            score = colors / omega
            if score > best_score + 1e-12 or (abs(score - best_score) < 1e-12 and len(seq) < len(best_seq)):
              best_seq = seq
              best_score = score
              best_meta = ("starts=%s depth=%d interleave=%s caps_before=%s caps=%s" % (
                              starts, depth, interleave, caps_before, caps), colors, omega, len(seq))

  if best_seq is None:
    best_seq = original()

  # Print the validation line
  if best_meta is None:
    best_meta = ("unknown", 0, 0, 0)
  best_desc, best_colors, best_omega, best_n = best_meta
  ratio = best_colors / best_omega if best_omega > 0 else 0.0
  print(f"construct_intervals: depth=0, order={best_desc}, n={best_n}, FirstFit={best_colors}, OPT={best_omega}, ratio={ratio:.3f}")

  return best_seq

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()