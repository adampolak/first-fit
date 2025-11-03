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
    """Return (colors_used, omega) for a given sequence tuple representation."""
    seq = list(seq_tuple)

    # Fast FirstFit coloring using per-color sorted-blocks with binary search
    color_blocks = []

    def can_place(block, new):
      l, r = new
      lo, hi = 0, len(block)
      while lo < hi:
        mid = (lo + hi) // 2
        if block[mid][0] < l:
          lo = mid + 1
        else:
          hi = mid
      pos = lo
      if pos > 0 and block[pos - 1][1] > l:
        return False
      if pos < len(block) and block[pos][0] < r:
        return False
      return True

    def insert_pos(block, l):
      lo, hi = 0, len(block)
      while lo < hi:
        mid = (lo + hi) // 2
        if block[mid][0] < l:
          lo = mid + 1
        else:
          hi = mid
      return lo

    for (l, r) in seq:
      placed = False
      for block in color_blocks:
        if can_place(block, (l, r)):
          pos = insert_pos(block, l)
          block.insert(pos, (l, r))
          placed = True
          break
      if not placed:
        color_blocks.append([(l, r)])
    colors_used = len(color_blocks)

    # Omega (clique number) via sweep-line on the full sequence
    events = []
    for (l, r) in seq:
      events.append((l, 1))
      events.append((r, -1))
    events.sort(key=lambda x: (x[0], x[1]))
    active = 0
    best = 0
    for _, delta in events:
      active += delta
      if active > best:
        best = active
    omega = best

    return (colors_used, omega)

  # -------------------- Orderings and adversarial machinery --------------------
  def order_left_first(seq):
    return sorted(seq, key=lambda x: (x[0], x[1]))

  def order_right_first(seq):
    return sorted(seq, key=lambda x: (-x[0], -x[1]))

  def order_short_first(seq):
    return sorted(seq, key=lambda x: ((x[1] - x[0]), x[0]))

  def order_long_first(seq):
    return sorted(seq, key=lambda x: (-(x[1] - x[0]), x[0]))

  def order_by_right_endpoint(seq):
    return sorted(seq, key=lambda x: (x[1], x[0]))

  def order_hybrid_blockwise(seq):
    # Partition into 4 blocks and apply different local orders
    n = len(seq)
    if n < 12:
      return list(seq)
    k = 4
    b = n // k
    parts = [list(seq[i * b:(i + 1) * b]) for i in range(k - 1)]
    parts.append(list(seq[(k - 1) * b:]))
    out = []
    out += parts[0]
    if len(parts) > 1:
      out += order_left_first(parts[1])
    if len(parts) > 2:
      out += list(reversed(parts[2]))
    if len(parts) > 3:
      out += order_short_first(parts[3])
    return out

  def adversarial_order(seq):
    """
    Adaptive adversary: iteratively pick an interval that blocks all current colors if possible,
    breaking ties by maximum coverage among unplaced intervals, then by shorter length.
    Uses bitsets; falls back for very large n.
    """
    n = len(seq)
    if n == 0:
      return []
    if n > 1200:
      # For very large instances, use a strong deterministic stress order
      return order_by_right_endpoint(seq)

    # Precompute intersection bitsets
    neighbors = [0] * n
    for i in range(n):
      li, ri = seq[i]
      for j in range(i + 1, n):
        lj, rj = seq[j]
        if li < rj and lj < ri:
          neighbors[i] |= (1 << j)
          neighbors[j] |= (1 << i)

    unplaced = set(range(n))
    color_sets = []   # list of bitsets, one per color
    assigned = [-1] * n
    ordering = []

    def assign_color(i):
      ni = neighbors[i]
      for c, cm in enumerate(color_sets):
        if (ni & cm) == 0:
          color_sets[c] |= (1 << i)
          assigned[i] = c
          return c
      # new color
      color_sets.append(1 << i)
      assigned[i] = len(color_sets) - 1
      return assigned[i]

    while unplaced:
      used = len(color_sets)
      # bitset for unplaced
      mask_unplaced = 0
      for j in unplaced:
        mask_unplaced |= (1 << j)
      best_i = None
      best_key = None
      for i in list(unplaced):
        ni = neighbors[i]
        # how many colors it intersects
        blocked = 0
        for cm in color_sets:
          if (ni & cm) != 0:
            blocked += 1
        new_color = 1 if blocked == used else 0
        coverage = ((ni & mask_unplaced) & ~(1 << i)).bit_count()
        length = seq[i][1] - seq[i][0]
        key = (new_color, blocked, coverage, -length)
        if best_key is None or key > best_key:
          best_key = key
          best_i = i
      ordering.append(best_i)
      assign_color(best_i)
      unplaced.remove(best_i)

    return [seq[i] for i in ordering]

  def best_ordered_seq(seq):
    """
    Try several deterministic orders and an adaptive adversarial order (budget-gated).
    Return (best_ordered_sequence, colors_used_by_FirstFit, omega).
    """
    orders = []
    orders.append(list(seq))  # as-built
    orders.append(list(reversed(seq)))
    orders.append(order_left_first(seq))
    orders.append(order_by_right_endpoint(seq))
    orders.append(order_short_first(seq))
    orders.append(order_long_first(seq))
    orders.append(order_hybrid_blockwise(seq))
    try:
      orders.append(adversarial_order(seq))
    except Exception:
      # safety net in case of pathologies
      pass
    best_seq = None
    best_colors = -1
    best_omega = 0
    for cand in orders:
      colors, omega = eval_seq_cached(tuple(cand))
      if omega <= 0:
        continue
      # primary objective: maximize colors/omega
      if best_omega == 0 or (colors / omega) > (best_colors / best_omega):
        best_seq = cand
        best_colors = colors
        best_omega = omega
    if best_seq is None:
      colors, omega = eval_seq_cached(tuple(seq))
      return list(seq), colors, omega
    return best_seq, best_colors, best_omega

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

  # Cap templates bridging adjacent and non-adjacent towers (deterministic scaffolds)
  def cap_templates_for_starts(starts):
    s = list(starts)
    t = []
    if len(s) >= 4:
      # Canonical pattern and a few balanced variants
      t.append([(s[0] - 1, s[0] + 3),
                (s[1] + 2, s[2] + 2),
                (s[2] - 1, s[3] + 1),
                (s[1] - 2, s[2] + 3)])
      t.append([(s[0], s[1] + 2),
                (s[1] + 1, s[2] + 3),
                (s[2], s[3] + 2),
                (s[0] + 2, s[2] + 1)])
      t.append([(s[0] - 1, s[1] + 3),
                (s[1] - 1, s[2] + 3),
                (s[2] - 1, s[3] + 3),
                (s[0] + 1, s[3] + 1)])
      t.append([(s[0] - 1, s[1] + 1),
                (s[2] - 1, s[3] + 1),
                (s[1], s[2] + 2),
                (s[0] + 2, s[3])])
    else:
      for i in range(len(s) - 1):
        t.append([(s[i], s[i + 1] + 2)] * 4)
    return t

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
    (4, 8, 12, 16),
    (2, 5, 9, 13),
    (2, 7, 12, 17),
    (3, 6, 9, 12),
    (2, 6, 10, 14, 18),
    (2, 6, 10, 14, 18, 22),
    (3, 7, 11, 15, 19),
    (4, 8, 12, 16, 20),
  ]
  depth_options = [3, 4, 5]
  trials_per_combo = 50

  # Baseline evaluation (with adversarial/best ordering)
  baseline = original()
  base_ordered, baseline_colors, baseline_omega = best_ordered_seq(baseline)
  if baseline_omega > 0:
    best_seq = base_ordered
    best_score = baseline_colors / baseline_omega
    best_meta = ("baseline+order", baseline_colors, baseline_omega, len(base_ordered))

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
  alt_ordered, alt_colors, alt_omega = best_ordered_seq(alt)
  if alt_omega > 0:
    alt_score = alt_colors / alt_omega
    if alt_score > best_score + 1e-12:
      best_seq = alt_ordered
      best_score = alt_score
      best_meta = ("alt_anchor=(3,7,11,15)+order", alt_colors, alt_omega, len(alt_ordered))

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
      # Deterministic template-based caps to couple towers
      templates = cap_templates_for_starts(starts)
      candidates += templates
      # Jittered variants around templates (budgeted)
      rng = random.Random(1337 + depth + len(starts))
      for tpl in templates:
        for _ in range(2):
          cc = []
          for (a, b) in tpl:
            da = rng.choice([-2, -1, 0, 1, 2])
            db = rng.choice([-1, 0, 1, 2, 3])
            aa = a + da
            bb = b + db
            if bb - aa < 3:
              bb = aa + 3
            cc.append((aa, bb))
          candidates.append(cc)
      # Random cap-sets to diversify search
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
            ordered_seq, colors, omega = best_ordered_seq(seq)
            if omega <= 0:
              continue
            score = colors / omega
            if score > best_score + 1e-12 or (abs(score - best_score) < 1e-12 and len(ordered_seq) < len(best_seq or [])):
              best_seq = ordered_seq
              best_score = score
              best_meta = ("starts=%s depth=%d interleave=%s caps_before=%s caps=%s" % (
                              starts, depth, interleave, caps_before, caps), colors, omega, len(ordered_seq))

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