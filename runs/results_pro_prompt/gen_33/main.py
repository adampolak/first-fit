# EVOLVE-BLOCK-START

def construct_intervals():
  """
  Construct a sequence of open intervals presented to FirstFit to maximize
  the ratio: FirstFit colors / clique-number (ω).
  This version blends fractal backbones, multiple seed/cap templates,
  a broad set of orderings, and adversarial-like refinement for better results.
  """

  # Utilities: FirstFit and clique number (ω) for open intervals
  def firstfit_colors(seq):
    end_times = []
    for l, r in seq:
      for i, et in enumerate(end_times):
        if et <= l:
          end_times[i] = r
          break
      else:
        end_times.append(r)
    return len(end_times)

  def clique_number(seq):
    events = []
    for l, r in seq:
      events.append((l, 1))
      events.append((r, -1))
    # open intervals: ends (-1) processed before starts (+1) at ties
    events.sort(key=lambda e: (e[0], e[1]))
    active = 0
    best = 0
    for _, delta in events:
      active += delta
      if active > best:
        best = active
    return best

  # -------------------- Ordering utilities --------------------

  def block_hybrid_order(seq):
    """Split seq into 4 blocks and reorder each block with a different heuristic."""
    n = len(seq)
    if n == 0:
      return []
    bs = (n + 3) // 4
    blocks = [seq[i*bs:(i+1)*bs] for i in range(4)]
    # pad if needed
    while len(blocks) < 4:
      blocks.append([])
    b0 = blocks[0]
    b1 = list(reversed(blocks[1]))
    b2 = sorted(blocks[2], key=lambda iv: (iv[1]-iv[0], iv[0]))
    b3 = sorted(blocks[3], key=lambda iv: (-(iv[1]-iv[0]), iv[0]))
    return b0 + b1 + b2 + b3

  def order_fractal_identity(seq):  # same as identity
    return list(seq)

  def order_reversed(seq):
    return list(reversed(seq))

  def order_left_first(seq):
    return sorted(seq, key=lambda x: (x[0], x[1]))

  def order_right_first(seq):
    return sorted(seq, key=lambda x: (-x[0], -x[1]))

  def order_short_first(seq):
    return sorted(seq, key=lambda x: (x[1] - x[0], x[0]))

  def order_long_first(seq):
    return sorted(seq, key=lambda x: (-(x[1] - x[0]), x[0]))

  def order_hybrid_blockwise(seq):
    """Block-wise reordering to diversify local structure."""
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

  def order_by_right_endpoint(seq):
    return sorted(seq, key=lambda x: (x[1], x[0]))

  def order_alt_short_long(seq):
    """Alternating small/large intervals to stress FirstFit."""
    n = len(seq)
    if n <= 2:
      return list(seq)
    by_len = sorted(seq, key=lambda x: (x[1] - x[0], x[0]))
    l, r = 0, n - 1
    take_left = True
    out = []
    while l <= r:
      if take_left:
        out.append(by_len[l]); l += 1
      else:
        out.append(by_len[r]); r -= 1
      take_left = not take_left
    return out

  ORDERERS = [
    ("fractal_identity", order_fractal_identity),
    ("reversed", order_reversed),
    ("left_first", order_left_first),
    ("right_first", order_right_first),
    ("short_first", order_short_first),
    ("long_first", order_long_first),
    ("hybrid_blockwise", order_hybrid_blockwise),
    ("right_endpoint", order_by_right_endpoint),
    ("alt_short_long", order_alt_short_long),
    ("hybrid_identity", block_hybrid_order),
  ]

  # -------------------- Adversarial-like ordering (deterministic) --------------------

  def adversarial_order(seq):
    """Deterministic greedy surrogate for an adversarial order: pick next
    interval with the largest degree (most overlaps with others)."""
    n = len(seq)
    if n == 0:
      return []
    # Build adjacency (overlap) graph
    overlaps = [set() for _ in range(n)]
    for i in range(n):
      li, ri = seq[i]
      for j in range(i + 1, n):
        lj, rj = seq[j]
        if li < rj and lj < ri:
          overlaps[i].add(j)
          overlaps[j].add(i)

    unplaced = set(range(n))
    order = []
    while unplaced:
      # pick a node with maximum degree (tie-breaker: larger interval length)
      best = max(unplaced, key=lambda i: (len(overlaps[i] & unplaced), (seq[i][1] - seq[i][0])))
      order.append(best)
      unplaced.remove(best)
    return [seq[i] for i in order]

  # -------------------- Fractal backbones and cap templates --------------------

  def build_fractal(depth, starts, caps, interleave=False, caps_before=False):
    """Recursive fractal-like backbone with given parameters."""
    T = [(0.0, 1.0)]
    for _ in range(depth):
      lo = min(l for l, _ in T)
      hi = max(r for _, r in T)
      delta = hi - lo
      S = []
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
    # Normalize
    lo = min(l for l, r in T)
    return [(l - lo, r - lo) for l, r in T]

  def cap_templates_for_starts(starts):
    """Produce several cap templates for a given starts tuple."""
    s = starts
    templates = []
    if len(s) >= 4:
      templates.append([
        (s[0] - 1, s[0] + 3),
        (s[1] + 2, s[2] + 2),
        (s[2] - 1, s[3] + 1),
        (s[1] - 2, s[2] + 3)
      ])
      templates.append([
        (s[0], s[1] + 2),
        (s[1] + 1, s[2] + 3),
        (s[2], s[3] + 2),
        (s[0] + 2, s[2] + 1)
      ])
      templates.append([
        (s[0] - 1, s[1] + 3),
        (s[1] - 1, s[2] + 3),
        (s[2] - 1, s[3] + 3),
        (s[0] + 1, s[3] + 1)
      ])
      templates.append([
        (s[0] - 1, s[1] + 1),
        (s[2] - 1, s[3] + 1),
        (s[1], s[2] + 2),
        (s[0] + 2, s[3])
      ])
    else:
      for i in range(len(s) - 1):
        templates.append([(s[i], s[i + 1] + 2)] * 4)
    return templates

  # -------------------- Candidate generation and screening --------------------

  class Candidate:
    __slots__ = ("intervals", "omega")
    def __init__(self, intervals, omega):
      self.intervals = intervals
      self.omega = omega

  def screen_candidates(cands, top_k=6):
    """Evaluate canonical orders and pick top_k by ratio colors/omega."""
    scored = []
    for cand in cands:
      if cand.omega <= 0:
        continue
      best_colors = -1
      best_seq = None
      for name, ord_fn in ORDERERS:
        seq_ord = ord_fn(list(cand.intervals))
        colors = firstfit_colors(seq_ord)
        if colors > best_colors:
          best_colors = colors
          best_seq = seq_ord
      ratio = best_colors / cand.omega
      scored.append((ratio, best_colors, cand.omega, len(cand.intervals), cand, best_seq))
    scored.sort(key=lambda x: (-x[0], x[3]))
    return scored[:top_k]

  def refine_with_adversary(top_entries):
    best_ratio = -1.0
    best_seq = None
    best_meta = None
    for ratio0, colors0, omega, n, cand, seq_ord in top_entries:
      # adversarial candidate
      seq_adv = adversarial_order(list(cand.intervals))
      colors_adv = firstfit_colors(seq_adv)
      ratio_adv = colors_adv / omega if omega > 0 else 0.0

      seq_adv2 = adversarial_order(list(seq_ord))
      colors_adv2 = firstfit_colors(seq_adv2)
      ratio_adv2 = colors_adv2 / omega if omega > 0 else 0.0

      candidates = [
        (ratio0, seq_ord, colors0),
        (ratio_adv, seq_adv, colors_adv),
        (ratio_adv2, seq_adv2, colors_adv2),
      ]
      for r, seq, cols in candidates:
        if r > best_ratio or (abs(r - best_ratio) < 1e-12 and (best_seq is None or len(seq) < len(best_seq))):
          best_ratio = r
          best_seq = seq
          best_meta = (cols, omega, len(seq))
    return best_seq, best_meta

  # -------------------- Parameter sweep (focused, deterministic) --------------------

  max_intervals = 800
  depths = [3, 4]
  starts_pool = [
    (2, 6, 10, 14),
    (3, 7, 11, 15),
    (4, 8, 12, 16),
    (2, 5, 9, 13),
    (3, 6, 9, 12),
  ]

  # Baseline seed
  baseline_caps = [(1, 5), (12, 16), (4, 9), (8, 13)]
  base_seq = build_fractal(4, (2, 6, 10, 14), baseline_caps, interleave=False, caps_before=False)
  base_omega = clique_number(base_seq)
  candidates = [Candidate(base_seq, base_omega)]

  # Seeded fractals
  seeds = [
    {'starts': (2, 6, 10, 14)},
    {'starts': (3, 7, 11, 15)},
    {'starts': (4, 8, 12, 16)},
  ]
  # For each seed, create several cap templates and depth variants
  for s in seeds:
    starts = s['starts']
    templates = cap_templates_for_starts(starts)
    # Build fractals for each template and depth
    for caps in templates:
      for depth in depths:
        seq = build_fractal(depth, starts, caps, interleave=False, caps_before=False)
        if 0 < len(seq) <= max_intervals:
          omega = clique_number(seq)
          if omega > 0:
            candidates.append(Candidate(seq, omega))

  # Explore with several starts patterns
  for starts in starts_pool:
    templates = cap_templates_for_starts(starts)
    for caps in templates:
      for depth in depths:
        for interleave in (False, True):
          for caps_before in (False, True):
            seq = build_fractal(depth, starts, caps, interleave=interleave, caps_before=caps_before)
            if 0 < len(seq) <= max_intervals:
              omega = clique_number(seq)
              if omega > 0:
                candidates.append(Candidate(seq, omega))

  # Screen and refine
  screened = screen_candidates(candidates, top_k=6)
  best_seq, best_meta = refine_with_adversary(screened)

  if best_seq is None:
    best_seq = base_seq

  # Final verification (optional, not printed)
  _ff = firstfit_colors(best_seq)
  _opt = clique_number(best_seq)
  return best_seq

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()