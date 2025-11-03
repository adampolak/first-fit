# EVOLVE-BLOCK-START

def construct_intervals(enable_alt_microphase=True):
  """
  Construct a sequence of intervals in FF order to maximize FirstFit/omega
  while keeping omega <= ~10 (target). Deterministic, CAP-aware pipeline.

  Args:
    enable_alt_microphase (bool): retain interface; enables the second micro-pass.

  Returns:
    List[(l, r)] integer open intervals (r > l).
  """

  # Global capacity guard (stay below evaluator hard bound)
  CAP = 9800
  BASE_SEED = 271828  # purely for deterministic choices; no RNG used

  # Six-template bank for the backbone rotation (recommendation 1)
  TEMPLATE_BANK_6 = [
    (2, 6, 10, 14),  # T0 classic KT
    (1, 5, 9, 13),   # T1 left-shifted
    (3, 7, 11, 15),  # T2 right-shifted
    (4, 8, 12, 16),  # T3 stretched-right
    (5, 9, 13, 17),  # T4 extra right-shift
    (2, 8, 12, 16),  # T5 inner symmetric
  ]

  # Stage helpers
  def span_delta(seq):
    lo = min(l for l, r in seq)
    hi = max(r for l, r in seq)
    d = hi - lo
    if d <= 0:
      d = 1
    return lo, hi, d

  def replicate_blocks(seq, starts):
    lo, hi, d = span_delta(seq)
    blocks = []
    for s in starts:
      base = s * d - lo
      blocks.append([(l + base, r + base) for (l, r) in seq])
    return blocks, d

  def interleave_blocks(blocks, order):
    out = []
    maxlen = max(len(b) for b in blocks) if blocks else 0
    for i in range(maxlen):
      for idx in order:
        blk = blocks[idx]
        if i < len(blk):
          out.append(blk[i])
    return out

  def append_kt_connectors(dst, starts, delta):
    s0, s1, s2, s3 = starts
    dst.append(((s0 - 1) * delta, (s1 - 1) * delta))  # left cap
    dst.append(((s2 + 2) * delta, (s3 + 2) * delta))  # right cap
    dst.append(((s0 + 2) * delta, (s2 - 1) * delta))  # cross 1
    dst.append(((s1 + 2) * delta, (s3 - 1) * delta))  # cross 2

  def normalize_and_cap(seq):
    if not seq:
      return []
    # Shift to nonnegative if needed and coerce to integers
    mn = min(l for l, r in seq)
    if mn < 0:
      seq = [(l - mn, r - mn) for (l, r) in seq]
    out = []
    for l, r in seq:
      li = int(l)
      ri = int(r)
      if ri <= li:
        ri = li + 1
      out.append((li, ri))
    if len(out) > CAP:
      out = out[:CAP]
    return out

  # Deterministic per-round spine densification (recommendation 5)
  def apply_density_multiplier(dst, center_frac=0.50, width_frac=0.20, mult=1.10, max_add=18):
    if not dst:
      return dst
    lo, hi, d = span_delta(dst)
    # Collect a central subset index range
    n = len(dst)
    c_idx = int(n * center_frac)
    halfw = max(1, int(n * width_frac * 0.5))
    left = max(0, c_idx - halfw)
    right = min(n, c_idx + halfw)
    subset = dst[left:right]

    # Shrink-and-duplicate a sparse sample from subset, capped by max_add and CAP
    if not subset:
      return dst
    stride = max(1, len(subset) // max(8, min(max_add, len(subset))))
    sample = subset[::stride][:max_add]
    eps = max(1, d // 1200)
    additions = []
    for (l, r) in sample:
      # Keep inside span and shrink slightly to avoid omega spikes
      nl = l + eps
      nr = r - eps
      if nr - nl > 0:
        additions.append((nl, nr))
    # Cap additions by remaining capacity
    room = CAP - len(dst)
    if room <= 0 or not additions:
      return dst
    if len(additions) > room:
      additions = additions[:room]
    dst.extend(additions)
    return dst

  # Two-pass micro-phase builder (recommendation 2)
  def build_micro(seq, budget, windows, add_long_cross=False, seed_boost=False):
    if not seq or budget <= 8:
      return []
    glo = min(l for l, r in seq)
    ghi = max(r for l, r in seq)
    G = max(1, ghi - glo)

    # Seed sizing depends on current length and pass type
    base_seed = max(8, min(36, len(seq) // 290))
    seed_sz = base_seed + (8 if seed_boost else 0)
    seed_sz = min(seed_sz, 48)

    stride = max(1, len(seq) // seed_sz)
    U = [seq[i] for i in range(0, len(seq), stride)][:seed_sz]
    if not U:
      return []
    ulo = min(l for l, r in U)

    # Build translated blocks for each window
    blocks = []
    for (fa, fb) in windows:
      fa = max(0.05, min(0.90, fa))
      fb = max(0.10, min(0.95, fb))
      win_lo = glo + int(round(fa * G))
      base = win_lo - ulo
      block = [(l + base, r + base) for (l, r) in U]
      # Alternate reversal by window parity to break symmetry
      tag = (int(round(fa * 100)) // 5) % 2
      if tag == 1:
        block = list(reversed(block))
      blocks.append(block)

    # Deterministic interleaving order based on base seed and window count
    order = list(range(len(blocks)))
    if (BASE_SEED + len(blocks)) % 3 == 0:
      order.reverse()
    micro = interleave_blocks(blocks, order)

    # Fractional-span connectors at micro scale
    micro_connectors = [
      (glo + int(round(0.08 * G)), glo + int(round(0.30 * G))),
      (glo + int(round(0.60 * G)), glo + int(round(0.92 * G))),
      (glo + int(round(0.26 * G)), glo + int(round(0.56 * G))),
      (glo + int(round(0.44 * G)), glo + int(round(0.78 * G))),
    ]
    for a, b in micro_connectors:
      if b > a:
        micro.append((a, b))

    # Optional single long cross (guarded)
    if add_long_cross:
      a = glo + int(round(0.18 * G))
      b = glo + int(round(0.84 * G))
      if b > a:
        micro.append((a, b))

    if len(micro) > budget:
      micro = micro[:budget]
    return micro

  # Post-spine deterministic long-range connectors (recommendation 3)
  def build_long_range_connectors(seq, at_most=6):
    if not seq or at_most <= 0:
      return []
    lo, hi, d = span_delta(seq)
    # Choose 4–6 connectors anchored to global span; endpoints reproducible from BASE_SEED
    picks = []
    anchors = [
      (0.12, 0.62),
      (0.22, 0.78),
      (0.38, 0.86),
      (0.48, 0.90),
      (0.30, 0.70),
      (0.15, 0.85),
    ]
    # Deterministic rotation
    rot = (BASE_SEED // 7) % len(anchors)
    for i in range(min(at_most, len(anchors))):
      a_frac, b_frac = anchors[(rot + i) % len(anchors)]
      L = lo + max(1, int(round(a_frac * d)))
      R = lo + max(1, int(round(b_frac * d)))
      if R > L:
        picks.append((L, R))
    # Guard return by CAP budget where appended
    return picks

  # Backbone: six rounds with deterministic interleaving schedule (recommendation 4)
  T = [(0, 1)]  # seed
  for ridx in range(6):
    # Predict size growth: sz -> 4*sz + 4
    nxt = 4 * len(T) + 4
    if nxt > CAP:
      break
    starts = TEMPLATE_BANK_6[(BASE_SEED + ridx) % len(TEMPLATE_BANK_6)]
    blocks, delta = replicate_blocks(T, starts)

    # Interleave on even rounds; sequential with reverse on odd rounds
    if ridx % 2 == 0:
      order = sorted(range(4), key=lambda i: starts[i])  # deterministic forward
      S = interleave_blocks(blocks, order)
    else:
      order = list(reversed(range(4)))  # reversed sequential
      S = []
      for i in order:
        S.extend(blocks[i])

    # Append KT connectors for this round
    append_kt_connectors(S, starts, delta)

    # Apply per-round density multiplier to a central subset (conservative)
    S = apply_density_multiplier(S, center_frac=0.50, width_frac=0.18, mult=1.10, max_add=14)

    # Capacity clamp after each round
    if len(S) > CAP:
      S = S[:CAP]

    T = S

  # Two-pass micro-phase budgeting anchored to the backbone (recommendation 2)
  # Pass-1 windows A; Pass-2 windows B, non-overlapping families
  windows_A = [(0.12, 0.22), (0.35, 0.45), (0.58, 0.68), (0.80, 0.90)]
  windows_B = [(0.06, 0.14), (0.28, 0.38), (0.54, 0.64), (0.74, 0.84)]

  # Budget split: allocate remaining CAP deterministically 55% for pass-1, 45% for pass-2
  room = CAP - len(T)
  if room > 8:
    bud1 = int(room * 0.55)
    bud2 = room - bud1
    # Build pass-1 (no long cross, moderate seed)
    micro1 = build_micro(T, bud1, windows_A, add_long_cross=False, seed_boost=False)
    if micro1:
      avail = CAP - len(T)
      if len(micro1) > avail:
        micro1 = micro1[:avail]
      T.extend(micro1)

  room = CAP - len(T)
  if enable_alt_microphase and room > 8:
    # Pass-2 (long cross enabled, slightly boosted seed)
    bud2 = room  # use remaining room
    micro2 = build_micro(T, bud2, windows_B, add_long_cross=True, seed_boost=True)
    if micro2:
      avail = CAP - len(T)
      if len(micro2) > avail:
        micro2 = micro2[:avail]
      T.extend(micro2)

  # Deterministic long-range connectors appended after micro-phases (recommendation 3)
  room = CAP - len(T)
  if room > 0:
    lr = build_long_range_connectors(T, at_most=min(6, room))
    if lr:
      add = lr[:room]
      T.extend(add)

  # Final normalization and capacity enforcement
  return normalize_and_cap(T)

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()