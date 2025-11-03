# EVOLVE-BLOCK-START

def construct_intervals():
  """
  Construct a sequence of open intervals presented to FirstFit.
  Enhanced search: broaden fractal seeds (starts/caps), adaptive depths,
  and diverse deterministic orderings including block-wise hybrids.

  Returns:
    intervals: list of tuples (l, r) representing open intervals (l, r)
  """
  import math
  import random
  from collections import defaultdict

  rng = random.Random(20251023)

  # ---------- Core utilities ----------

  def firstfit_for_order(intervals):
    # Standard FirstFit in the given order (open intervals)
    colors = []
    for i, (l, r) in enumerate(intervals):
      used = set()
      for j, (lj, rj) in enumerate(intervals[:i]):
        if lj < r and l < rj:
          used.add(colors[j])
      c = 1
      while c in used:
        c += 1
      colors.append(c)
    return colors

  def firstfit_colors(intervals):
    assigned = firstfit_for_order(intervals)
    return max(assigned) if assigned else 0

  def clique_number(intervals):
    # Sweep-line for open intervals: end events before start events on ties
    events = []
    for l, r in intervals:
      events.append((l, 1))   # start
      events.append((r, -1))  # end
    events.sort(key=lambda x: (x[0], x[1]))
    curr = 0
    best = 0
    for _, delta in events:
      curr += delta
      if curr > best:
        best = curr
    return best

  def normalize(seq):
    if not seq:
      return seq
    lo = min(l for l, r in seq)
    return [(l - lo, r - lo) for (l, r) in seq]

  # ---------- Fractal builder with meta tags ----------

  # Each item is (l, r, tag, root), where
  #  - tag in {'rep','cap','base'}
  #  - root is the top-level replicate index (for block-wise hybrids)
  def build_fractal_labelled(depth, starts, caps, order_caps_after=True):
    # Initial unit tile
    T = [(0.0, 1.0, 'base', -1, 0)]  # (l, r, tag, root, level)
    for level in range(1, depth + 1):
      lo = min(l for l, r, _, __, ___ in T)
      hi = max(r for l, r, _, __, ___ in T)
      delta = hi - lo
      S = []
      # replicate blocks
      for idx, start in enumerate(starts):
        shift = delta * start - lo
        for (l, r, tag, root, lv) in T:
          # root is set to current top-level replicate index at first expansion
          new_root = idx if level == 1 else root
          S.append((shift + l, shift + r, 'rep' if tag != 'cap' else 'cap', new_root, level))
      # add caps at this level
      for (a, b) in caps:
        S.append((delta * a, delta * b, 'cap', -1, level))
      T = S

    # normalize
    base = min(l for l, r, _, __, ___ in T)
    T = [(l - base, r - base, tag, root, lv) for (l, r, tag, root, lv) in T]

    # Provide a default "fractal" order: as built
    return T

  # ---------- Cap generators for starts ----------

  def caps_adjacent(starts, pad_left=0.0, pad_right=1.0, include_extremes=True):
    # Bridge between consecutive replicate windows [start, start+1]
    s = list(starts)
    s.sort()
    caps = []
    for i in range(len(s) - 1):
      caps.append((s[i] + pad_left, s[i + 1] + pad_right))
    if include_extremes:
      # extra l/r extremes
      caps.append((s[0] - 1.0, s[0] + 3.0))
      caps.append((s[-1] - 3.0, s[-1] + 1.0))
    return caps

  def caps_sliding(starts, window=2, pad=(0.0, 1.0)):
    # Sliding windows of size=window over starts
    s = list(starts)
    s.sort()
    caps = []
    for i in range(0, len(s) - window + 1):
      L = s[i] + pad[0]
      R = s[i + window - 1] + 1.0 + pad[1]
      caps.append((L, R))
    return caps

  def generate_anchor_sets():
    # Multiple starts with different steps and offsets; include original
    starts_list = [
      (2, 6, 10, 14),         # original seed (M=4)
      (3, 7, 11, 15),         # shift +1
      (4, 8, 12, 16),         # shift +2
      (2, 5, 8, 11, 14),      # M=5, step 3
      (2, 6, 10, 14, 18),     # M=5, step 4
      (2, 6, 10, 14, 18, 22), # M=6, step 4
    ]
    seeds = []

    for starts in starts_list:
      M = len(starts)
      # Family A: adjacent caps with common paddings
      for (pl, pr) in [(0.0, 1.0), (1.0, 1.0), (-0.5, 1.5), (0.0, 2.0)]:
        capsA = caps_adjacent(starts, pl, pr, include_extremes=True)
        seeds.append(("adj", starts, capsA))
      # Family B: sliding caps of size 2 and 3 where possible
      for window in [2, 3]:
        for pad in [(0.0, 1.0), (0.5, 1.0), (0.0, 1.5)]:
          capsB = caps_sliding(starts, window=window, pad=pad)
          # Lightly augment with one extreme when M small
          if M <= 4 and capsB:
            capsB = list(capsB) + [(starts[0] - 1.0, starts[0] + 3.0)]
          seeds.append((f"slide{window}", starts, capsB))
      # Family C: original paper-like caps when M=4 with symmetric pattern
      if starts == (2, 6, 10, 14):
        seeds.append(("paper", starts, [(1, 5), (12, 16), (4, 9), (8, 13)]))

    # Add a few randomized variants around each deterministic seed to enrich diversity
    enriched = []
    for name, starts, caps in seeds:
      enriched.append((name, starts, caps))
      # jitter caps slightly in a deterministic way
      for j in range(min(2, max(0, 6 - len(caps)))):  # up to 2 jitters
        jittered = []
        local_rng = random.Random((len(starts), len(caps), j, 17))
        for (a, b) in caps:
          da = (local_rng.random() - 0.5) * 0.4  # +/- 0.2
          db = (local_rng.random() - 0.5) * 0.4
          if b + db > a + da + 0.3:  # keep open interval valid and not too tiny
            jittered.append((a + da, b + db))
          else:
            jittered.append((a, b))
        enriched.append((name + "_jit", starts, jittered))
    return enriched

  # ---------- Orders ----------

  def apply_order_labelled(labelled, order_name):
    # labelled: list of (l, r, tag, root, level)
    if order_name == "fractal":
      return labelled[:]  # construction order
    if order_name == "reversed":
      return list(reversed(labelled))
    if order_name == "short_first":
      return sorted(labelled, key=lambda x: ((x[1] - x[0]), x[0]))
    if order_name == "long_first":
      return sorted(labelled, key=lambda x: (-(x[1] - x[0]), x[0]))
    if order_name == "left_first":
      return sorted(labelled, key=lambda x: (x[0], x[1]))
    if order_name == "right_first":
      return sorted(labelled, key=lambda x: (-x[0], x[1]))
    if order_name == "cap_first":
      caps = [t for t in labelled if t[2] == 'cap']
      reps = [t for t in labelled if t[2] != 'cap']
      return caps + reps
    if order_name == "cap_last":
      caps = [t for t in labelled if t[2] == 'cap']
      reps = [t for t in labelled if t[2] != 'cap']
      return reps + caps
    if order_name.startswith("block_hybrid"):
      # Partition by top-level root id into up to 4 blocks, apply different sub-orders
      blocks = defaultdict(list)
      caps = []
      for item in labelled:
        if item[2] == 'cap':
          caps.append(item)
        else:
          blocks[item[3]].append(item)  # group by root
      # sort block keys
      keys = sorted(blocks.keys())
      # split into 4 slices
      chunks = [keys[i::4] for i in range(4)]
      B1 = [x for k in chunks[0] for x in sorted(blocks[k], key=lambda y: (y[0], y[1]))]  # left_first
      B2 = [x for k in chunks[1] for x in sorted(blocks[k], key=lambda y: (-(y[1] - y[0]), y[0]))]  # long_first
      B3 = [x for k in chunks[2] for x in sorted(blocks[k], key=lambda y: (-y[0], y[1]))]  # right_first
      B4 = [x for k in chunks[3] for x in sorted(blocks[k], key=lambda y: ((y[1] - y[0]), y[0]))]  # short_first
      if order_name.endswith("_caps_first"):
        return caps + B1 + B2 + B3 + B4
      elif order_name.endswith("_caps_last"):
        return B1 + B2 + B3 + B4 + caps
      else:
        # interleave caps between blocks
        return B1 + caps + B2 + B3 + B4
    return labelled[:]  # default

  def strip_pairs(labelled):
    return [(l, r) for (l, r, _, __, ___) in labelled]

  # ---------- Depth budget ----------

  def max_depth_for_budget(M, K, budget):
    # N(d) = M^d + K * (M^d - 1) / (M - 1)
    # Find largest d with N(d) <= budget (M>=2)
    d = 0
    while True:
      Md = M ** d
      N = Md + (K * (Md - 1) // (M - 1) if M > 1 else Md + K * d)
      if N > budget:
        return max(0, d - 1)
      d += 1

  # ---------- Candidate search ----------

  seeds = generate_anchor_sets()

  orders_to_try = [
    "fractal",
    "reversed",
    "short_first",
    "long_first",
    "left_first",
    "right_first",
    "cap_first",
    "cap_last",
    "block_hybrid",
    "block_hybrid_caps_first",
    "block_hybrid_caps_last",
  ]

  budget_N = 800
  best_seq_pairs = None
  best_meta = None
  best_ratio = -1.0

  # Limit the total number of evaluated candidates
  max_evals = 220
  eval_count = 0

  # Deterministically shuffle seeds to avoid always starting with the same
  rng.shuffle(seeds)

  for seed_idx, (seed_name, starts, caps) in enumerate(seeds):
    M = len(starts)
    K = len(caps)
    # limit depths to practical range and budget
    dmax = max_depth_for_budget(M, K, budget_N)
    # We want depth >= 2; if too small, skip
    depths = []
    for d in range(2, min(5, dmax) + 1):
      depths.append(d)
    if not depths:
      continue

    # Try a couple of depths (favor deeper first)
    for depth in sorted(depths, reverse=True):
      labelled = build_fractal_labelled(depth, starts, caps)
      n = len(labelled)
      if n > budget_N:
        continue

      # precompute base opt once for this set (independent of order)
      seq_pairs_full = strip_pairs(labelled)
      opt_val = max(1, clique_number(seq_pairs_full))  # avoid div-by-zero

      # Orders loop
      for order_name in orders_to_try:
        candidate_labelled = apply_order_labelled(labelled, order_name)
        candidate = strip_pairs(candidate_labelled)

        # Evaluate FirstFit
        ff_colors = firstfit_colors(candidate)
        ratio = ff_colors / opt_val

        # Track best; tie-breaker: higher ff_colors, then fewer intervals
        improved = False
        if ratio > best_ratio + 1e-12:
          improved = True
        elif abs(ratio - best_ratio) <= 1e-12:
          # tie-breaker policies
          prev_colors = best_meta[2] if best_meta else -1
          prev_n = best_meta[5] if best_meta else 10**9
          if ff_colors > prev_colors or (ff_colors == prev_colors and n < prev_n):
            improved = True

        if improved:
          best_ratio = ratio
          best_seq_pairs = candidate
          best_meta = (depth, order_name, ff_colors, opt_val, seed_name, n, M, K, starts, caps)

        eval_count += 1
        if eval_count >= max_evals:
          break
      if eval_count >= max_evals:
        break
    if eval_count >= max_evals:
      break

  # Fallback to original if nothing found (should not happen)
  if best_seq_pairs is None:
    def old_build():
      T = [(0.0, 1.0)]
      for _ in range(2):
        lo = min(l for l, r in T)
        hi = max(r for l, r in T)
        delta = hi - lo
        S = []
        for start in (2, 6, 10, 14):
          S += [(delta * start + l - lo, delta * start + r - lo) for l, r in T]
        S += [(delta * 1, delta * 5), (delta * 12, delta * 16), (delta * 4, delta * 9), (delta * 8, delta * 13)]
        T = S
      return T
    best_seq_pairs = old_build()
    opt_val = clique_number(best_seq_pairs)
    ff_colors = firstfit_colors(best_seq_pairs)
    best_meta = (2, "fallback", ff_colors, opt_val, "paper", len(best_seq_pairs), 4, 4, (2,6,10,14), [(1,5),(12,16),(4,9),(8,13)])
    best_ratio = ff_colors / max(1, opt_val)

  # Verification and normalization
  seq = normalize(best_seq_pairs)
  verified_ff = firstfit_colors(seq)
  verified_opt = max(1, clique_number(seq))
  # Short validation print
  depth, order_name, ff_colors, opt_val, seed_name, n_intervals, M, K, starts, caps = best_meta
  print(f"construct_intervals: depth={depth}, order={order_name}, seed={seed_name}, M={M}, K={K}, n={n_intervals}, FF={ff_colors}, OPT={opt_val}, ratio={best_ratio:.3f}")
  print(f"verification: FF={verified_ff}, OPT={verified_opt}, ratio={verified_ff/verified_opt:.3f}")
  return seq

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()