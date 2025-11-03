# EVOLVE-BLOCK-START

def construct_intervals():
  """
  Improved fractal + adversary pipeline that returns an arrival-ordered list
  of open intervals intended to maximize FirstFit / OPT (clique) ratio.
  """
  import random
  from functools import lru_cache
  import math
  random.seed(0)

  # -------------------- Basic utilities --------------------

  Interval = tuple  # (l, r)

  def sweep_clique(intervals):
    """Sweep-line for open intervals: process ends before starts at equal coords."""
    events = []
    for (l, r) in intervals:
      events.append((l, 1))
      events.append((r, -1))
    events.sort(key=lambda x: (x[0], x[1]))
    curr = best = 0
    for _, d in events:
      curr += d
      if curr > best:
        best = curr
    return best

  # -------------------- Fractal / tile builder --------------------

  def build_fractal(depth, starts, caps, interleave=False, caps_before=False):
    """
    Recursive fractal-like generator producing a set of intervals.
    - starts: tuple of integer offsets for block copies per level
    - caps: list of (a,b) intervals appended each level (scaled by delta)
    - interleave/caps_before control copy and cap ordering
    """
    T = [(0.0, 1.0)]
    for _ in range(depth):
      lo = min(l for l, r in T)
      hi = max(r for l, r in T)
      delta = hi - lo
      S = []
      if caps_before:
        for a, b in caps:
          S.append((delta * a, delta * b))
      if interleave:
        # round-robin interleaving of blocks (increases cross-block interference)
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
    # Normalize leftmost to 0 for numerical stability
    lo = min(l for l, r in T) if T else 0.0
    return [(l - lo, r - lo) for l, r in T]

  # -------------------- Cap templates --------------------

  def cap_templates_for_starts(starts):
    s = list(starts)
    tpls = []
    m = len(s)
    if m >= 4:
      tpls.append([(s[0] - 1, s[0] + 3),
                   (s[1] + 2, s[2] + 2),
                   (s[2] - 1, s[3] + 1),
                   (s[1] - 2, s[2] + 3)])
      tpls.append([(s[0], s[1] + 2),
                   (s[1] + 1, s[2] + 3),
                   (s[2], s[3] + 2),
                   (s[0] + 2, s[2] + 1)])
      tpls.append([(s[0] - 1, s[1] + 3),
                   (s[1] - 1, s[2] + 3),
                   (s[2] - 1, s[3] + 3),
                   (s[0] + 1, s[3] + 1)])
      tpls.append([(s[0] - 1, s[1] + 1),
                   (s[2] - 1, s[3] + 1),
                   (s[1], s[2] + 2),
                   (s[0] + 2, s[3])])
    else:
      for i in range(max(1, m - 1)):
        tpls.append([(s[i], s[min(i + 1, m - 1)] + 2)] * 4)
    return tpls

  # -------------------- Overlap bitsets & FirstFit simulation --------------------

  def compute_neighbors(seq):
    """Return list of bitset neighbors: neighbors[i] has bit j set iff intervals i and j overlap."""
    n = len(seq)
    neigh = [0] * n
    for i in range(n):
      li, ri = seq[i]
      # j > i loop to avoid double-check
      for j in range(i + 1, n):
        lj, rj = seq[j]
        if li < rj and lj < ri:
          neigh[i] |= 1 << j
          neigh[j] |= 1 << i
    return neigh

  def firstfit_count_with_neighbors(neigh, order_indices):
    """Simulate FirstFit using precomputed neighbor bitsets on a given order of indices."""
    color_masks = []  # list of ints (bitsets)
    for idx in order_indices:
      ni = neigh[idx]
      placed = False
      # find the smallest color that has no conflict
      for c in range(len(color_masks)):
        if (ni & color_masks[c]) == 0:
          color_masks[c] |= (1 << idx)
          placed = True
          break
      if not placed:
        color_masks.append(1 << idx)
    return len(color_masks)

  def firstfit_count_seq(seq):
    """Full FirstFit on seq in given order (seq is list of intervals); returns color count."""
    # compute neighbors for this sequence's indexing
    neigh = compute_neighbors(seq)
    order_indices = list(range(len(seq)))
    return firstfit_count_with_neighbors(neigh, order_indices)

  # -------------------- Strong adversarial ordering --------------------

  def adversarial_order_max_index(seq, sample_limit=250, rnd_seed=1234):
    """
    Greedy adversary: at each step select the interval that (when assigned now)
    would get the largest color index (i.e., prefers choices that create a new color).
    Tie-break by coverage with remaining intervals and by interval length (shorter preferred).
    Uses bitsets for speed.
    """
    n = len(seq)
    if n == 0:
      return []

    neigh = compute_neighbors(seq)
    full_mask = (1 << n) - 1
    unpicked_mask = full_mask
    color_masks = []  # masks of indices assigned to each color
    order_idx = []
    rng = random.Random(rnd_seed + n)

    # Precomputed lengths
    lengths = [seq[i][1] - seq[i][0] for i in range(n)]

    while unpicked_mask:
      remaining = unpicked_mask.bit_count()
      # candidates sampling
      if remaining <= sample_limit:
        picks = []
        m = unpicked_mask
        while m:
          lsb = (m & -m)
          i = (lsb.bit_length() - 1)
          picks.append(i)
          m &= m - 1
      else:
        # sample random indices from unpicked_mask
        picks = []
        # produce list of unpicked indices for sampling
        up_list = []
        m = unpicked_mask
        while m:
          lsb = (m & -m)
          up_list.append(lsb.bit_length() - 1)
          m &= m - 1
        picks = rng.sample(up_list, min(sample_limit, len(up_list)))

      best_i = None
      best_key = None

      mask_unpicked = unpicked_mask
      for i in picks:
        ni = neigh[i]
        # find assigned index if placed now
        assigned_idx = None
        for c, cmask in enumerate(color_masks):
          if (ni & cmask) == 0:
            assigned_idx = c
            break
        if assigned_idx is None:
          assigned_idx = len(color_masks)  # would create a new color

        # coverage: how many remaining intervals would it touch
        cov = (ni & (mask_unpicked & ~(1 << i))).bit_count()
        # blocked color count
        blocked = 0
        for cmask in color_masks:
          if (ni & cmask) != 0:
            blocked += 1
        # Prefer higher assigned_idx, then higher coverage, then more blocked, then shorter length
        key = (assigned_idx, cov, blocked, -lengths[i])
        if best_key is None or key > best_key:
          best_key = key
          best_i = i

      # Place best_i
      bi = best_i
      if bi is None:
        # fallback: take any unpicked bit
        bi = (unpicked_mask & -unpicked_mask).bit_length() - 1

      # assign to first-fit color index according to current color_masks
      ni = neigh[bi]
      placed = False
      for c in range(len(color_masks)):
        if (ni & color_masks[c]) == 0:
          color_masks[c] |= (1 << bi)
          placed = True
          break
      if not placed:
        color_masks.append(1 << bi)

      order_idx.append(bi)
      unpicked_mask &= ~(1 << bi)

    # translate indices into intervals
    return [seq[i] for i in order_idx]

  # -------------------- Candidate generation/search --------------------

  # Parameter guards
  max_intervals = 800
  depths = [3, 4, 5]
  # anchor families: include merged families to diversify fractal structure
  starts_families = [
    (2, 6, 10, 14),
    (3, 7, 11, 15),
    (4, 8, 12, 16),
    (2, 3, 6, 7, 10, 11, 14, 15),  # merged anchors
    (2, 4, 6, 8, 10, 12, 14, 16),
  ]
  # baseline caps
  baseline_caps = [(1, 5), (12, 16), (4, 9), (8, 13)]

  candidates = []
  # Always include the original baseline depth-4 variant as a safe candidate
  baseline = build_fractal(4, (2, 6, 10, 14), baseline_caps, interleave=False, caps_before=False)
  if 0 < len(baseline) <= max_intervals:
    candidates.append(("baseline", baseline))

  rng = random.Random(1337)
  candidate_limit = 150

  for starts in starts_families:
    templates = cap_templates_for_starts(starts)
    # generate jittered variants around templates (bounded number)
    caps_sets = []
    for tpl in templates:
      caps_sets.append(tpl)
      for _ in range(2):
        jitter = []
        for (a, b) in tpl:
          da = rng.choice([-2, -1, 0, 1, 2])
          db = rng.choice([-1, 0, 1, 2, 3])
          aa, bb = a + da, b + db
          if bb - aa < 3:
            bb = aa + 3
          jitter.append((aa, bb))
        caps_sets.append(jitter)

    # also include the baseline caps
    caps_sets.append(baseline_caps)

    for depth in depths:
      # rough size estimate: s * prev + p each level
      for caps in caps_sets:
        # estimate explosion and prune
        est = 1
        p_est = max(4, len(caps))
        feasible = True
        for _ in range(depth):
          est = len(starts) * est + p_est
          if est > max_intervals:
            feasible = False
            break
        if not feasible:
          continue
        for interleave in (False, True):
          for caps_before in (False, True):
            seq = build_fractal(depth, starts, caps, interleave=interleave, caps_before=caps_before)
            if 0 < len(seq) <= max_intervals:
              key = (len(seq), round(min(l for l, r in seq), 6), round(max(r for l, r in seq), 6))
              candidates.append((f"starts={starts} depth={depth} inter={interleave} caps_before={caps_before} caps={caps[:2]}", seq))
              if len(candidates) >= candidate_limit:
                break
        if len(candidates) >= candidate_limit:
          break
      if len(candidates) >= candidate_limit:
        break
    if len(candidates) >= candidate_limit:
      break

  # De-duplicate by (size, leftmost, rightmost)
  unique = {}
  for name, seq in candidates:
    if not seq:
      continue
    key = (len(seq), round(min(l for l, r in seq), 6), round(max(r for l, r in seq), 6))
    if key not in unique:
      unique[key] = (name, seq)
  candidates = list(unique.values())

  # -------------------- Evaluation & selection --------------------

  best_seq = None
  best_ratio = -1.0
  best_meta = None

  # small set of deterministic orderers to test as sanity
  def order_identity(x): return list(x)
  def order_reversed(x): return list(reversed(x))
  def order_left_first(x): return sorted(x, key=lambda iv: (iv[0], iv[1]))
  def order_short_first(x): return sorted(x, key=lambda iv: (iv[1] - iv[0], iv[0]))
  orderers = [("identity", order_identity), ("reversed", order_reversed), ("left", order_left_first), ("short", order_short_first)]

  # Evaluate each candidate using the strong adversary (plus a few deterministic orders)
  for idx, (name, seq) in enumerate(candidates):
    if not seq:
      continue
    omega = sweep_clique(seq)
    if omega <= 0:
      continue

    # 1) adversarial order (primary)
    try:
      seq_adv = adversarial_order_max_index(seq, sample_limit=250, rnd_seed=1000 + idx)
    except Exception:
      seq_adv = seq[:]
    # compute firstfit on seq_adv using neighbors of seq_adv
    ff_adv = firstfit_count_seq(seq_adv)
    ratio_adv = ff_adv / omega

    # 2) quick deterministic screening (best among a few)
    best_det_colors = -1
    best_det_seq = None
    for oname, ofn in orderers:
      ord_seq = ofn(seq)
      c = firstfit_count_seq(ord_seq)
      if c > best_det_colors:
        best_det_colors = c
        best_det_seq = ord_seq
    ratio_det = best_det_colors / omega

    # 3) also evaluate adversary when applied to the best deterministic sequence as seed
    try:
      seq_adv2 = adversarial_order_max_index(best_det_seq, sample_limit=200, rnd_seed=2000 + idx)
    except Exception:
      seq_adv2 = best_det_seq[:]
    ff_adv2 = firstfit_count_seq(seq_adv2)
    ratio_adv2 = ff_adv2 / omega

    # choose the best ordering among evaluated
    candidate_list = [
      (ratio_adv, ff_adv, seq_adv),
      (ratio_det, best_det_colors, best_det_seq),
      (ratio_adv2, ff_adv2, seq_adv2),
    ]
    candidate_list.sort(key=lambda x: (-x[0], -x[1], len(x[2])))

    cand_ratio, cand_colors, cand_seq = candidate_list[0][0], candidate_list[0][1], candidate_list[0][2]
    # Prefer higher ratio; tiebreaker: fewer intervals (smaller witness)
    if cand_ratio > best_ratio + 1e-12 or (abs(cand_ratio - best_ratio) < 1e-12 and len(cand_seq) < (len(best_seq) if best_seq else 10**9)):
      best_ratio = cand_ratio
      best_seq = cand_seq
      best_meta = (name, len(seq), cand_colors, omega)

  # Fallback to baseline if nothing found
  if best_seq is None:
    best_seq = baseline
    best_meta = ("fallback", len(baseline), firstfit_count_seq(baseline), sweep_clique(baseline))

  # Final cross-verification (sweep-line ω)
  desc, n_before, ff_before, omega_before = best_meta
  omega_check = sweep_clique(best_seq)
  ff_check = firstfit_count_seq(best_seq)
  ratio_check = ff_check / (omega_check if omega_check > 0 else 1)
  # Print a concise validation line (helpful during debugging/tuning)
  print(f"construct_intervals: candidate={desc}, n={len(best_seq)}, FirstFit={ff_check}, OPT={omega_check}, ratio={ratio_check:.3f}")

  return best_seq

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()