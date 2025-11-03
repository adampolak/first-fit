# EVOLVE-BLOCK-START

def construct_intervals(rounds=6,
                        rotate_starts=True,
                        reverse_block_parity=True,
                        interleave_blocks=True,
                        phase2_iters=1):
  """
  Deterministic KT-style spine with a rotating-template scaffold,
  enhanced by dual micro-waves and an adaptive coverage-aware corridor phase.

  Returns:
    intervals: list of (l, r) integer tuples, open intervals, in FF presentation order.
  """

  # Hard capacity guard to keep the total count < 10000
  CAP = 9800

  # Six strong start-pattern templates (diversified KT starts).
  template_bank6 = [
    (2, 6, 10, 14),  # classic KT
    (1, 5, 9, 13),   # left-shifted
    (3, 7, 11, 15),  # right-shifted
    (4, 8, 12, 16),  # stretched-right
    (2, 4, 8, 12),   # compressed-left pair
    (2, 8, 10, 12),  # inner-symmetric
  ]

  # Seed with one unit interval to allow six KT rounds within CAP.
  T = [(0, 1)]

  # Predictive size accounting to cap the number of full rounds (KT growth: size -> 4*size + 4).
  def round_next_size(sz):
    return 4 * sz + 4

  def max_rounds_within_cap(initial_size, max_rounds):
    sz = initial_size
    done = 0
    for _ in range(max(0, int(max_rounds))):
      nxt = round_next_size(sz)
      if nxt > CAP:
        break
      sz = nxt
      done += 1
    return done, sz

  def _span_delta(current_T):
    lo = min(l for l, r in current_T)
    hi = max(r for l, r in current_T)
    delta = hi - lo
    if delta <= 0:
      delta = 1
    return lo, hi, delta

  def _append_connectors(S, starts, delta, extra_cross=False):
    # Classic four connectors; preserves strong FF pressure while keeping omega modest.
    s0, s1, s2, s3 = starts
    S.append(((s0 - 1) * delta, (s1 - 1) * delta))  # left cap
    S.append(((s2 + 2) * delta, (s3 + 2) * delta))  # right cap
    S.append(((s0 + 2) * delta, (s2 - 1) * delta))  # cross 1
    S.append(((s1 + 2) * delta, (s3 - 1) * delta))  # cross 2
    # Optional long-range connector (used sparingly to avoid omega spikes)
    if extra_cross:
      S.append(((s0 + 3) * delta, (s3 + 3) * delta))

  def apply_round(current_T, starts, do_interleave=False, reverse_order=False, snake=False, extra_cross=False):
    lo, _, delta = _span_delta(current_T)

    # Build four translated blocks (snake reversal alternates internal order)
    blocks = []
    for idx, s in enumerate(starts):
      base = s * delta - lo
      block = [(l + base, r + base) for (l, r) in current_T]
      if snake and (idx % 2 == 1):
        block = list(reversed(block))
      blocks.append(block)

    # Build S either interleaving or sequential
    S = []
    if do_interleave:
      maxlen = max(len(b) for b in blocks)
      order = list(range(4))
      if reverse_order:
        order.reverse()
      for i in range(maxlen):
        for j in order:
          blk = blocks[j]
          if i < len(blk):
            S.append(blk[i])
    else:
      if reverse_order:
        blocks = list(reversed(blocks))
      for blk in blocks:
        S.extend(blk)

    # Append connectors (with optional extra long cross)
    _append_connectors(S, starts, delta, extra_cross=extra_cross)
    return S

  def _select_template(ridx):
    if not rotate_starts:
      return template_bank6[0]
    # Deterministic, diverse selection over six templates
    idx = (1729 * (ridx + 1) + 2654435761) % len(template_bank6)
    return template_bank6[idx]

  # Stage 1: KT spine with rotating templates, parity-based interleaving, and snake reversal.
  depth, _ = max_rounds_within_cap(len(T), rounds)
  for ridx in range(depth):
    starts = _select_template(ridx)
    do_inter = bool(interleave_blocks and (ridx % 2 == 0))   # interleave on even rounds
    rev = bool(reverse_block_parity and (ridx % 2 == 1))     # reverse block order on odd rounds
    snake = (ridx % 3 == 1)                                  # snake internal reversal on 1,4,...
    extra_cross = (ridx % 3 == 2)                            # rare long cross on 2,5,...
    T = apply_round(T, starts, do_interleave=do_inter, reverse_order=rev, snake=snake, extra_cross=extra_cross)

  # If we are close to capacity, return the strong baseline.
  if len(T) >= CAP - 8:
    return T

  # Utility: insert intervals near the tail (heightens late FF blocking)
  def _insert_near_tail(seq, intervals):
    out = list(seq)
    for i, iv in enumerate(intervals):
      pos = len(out) - (i * 2 + 1)
      if pos < 0:
        out.append(iv)
      else:
        out.insert(pos, iv)
    return out

  # Edge anchors to couple colors near boundaries without inflating omega.
  lo, hi, delta = _span_delta(T)
  edge_caps = []
  L1 = lo + max(1, int(round(0.05 * delta)))
  R1 = lo + max(2, int(round(0.20 * delta)))
  if R1 > L1:
    edge_caps.append((L1, R1))
  L2 = lo + max(1, int(round(0.80 * delta)))
  R2 = lo + max(2, int(round(0.95 * delta)))
  if R2 > L2:
    edge_caps.append((L2, R2))
  room = CAP - len(T)
  if room > 0 and edge_caps:
    T = _insert_near_tail(T, edge_caps[:room])

  # Recompute span after edge anchors and place spine caps.
  lo, hi, delta = _span_delta(T)

  def cap_at(a_frac, b_frac):
    L = lo + max(1, int(round(a_frac * delta)))
    R = lo + max(1, int(round(b_frac * delta)))
    if R <= L:
      R = L + 1
    return (L, R)

  # Long caps positioned as fractions of the current span; ensure monotone endpoints.
  caps = [cap_at(0.08, 0.60), cap_at(0.25, 0.75), cap_at(0.75, 0.92)]
  room = CAP - len(T)
  if room > 0:
    T = _insert_near_tail(T, caps[:room])

  # Additional interior spine caps (pairwise disjoint windows) to increase mixing safely
  lo, hi, delta = _span_delta(T)
  span = max(1, hi - lo)
  def frac_iv(a, b):
    L = lo + max(1, int(round(a * span)))
    R = lo + max(1, int(round(b * span)))
    return (L, R) if R > L else (L, L + 1)

  spine_caps = [
    frac_iv(0.06, 0.18),
    frac_iv(0.28, 0.40),
    frac_iv(0.46, 0.58),
    frac_iv(0.64, 0.76),
  ]
  room = CAP - len(T)
  if room > 0 and spine_caps:
    # Append (not near-tail) to diversify interaction times
    add = spine_caps[:room]
    T.extend(add)

  if len(T) >= CAP - 16:
    return T

  # Stage 2: delta2-driven micro waves using thin evenly-spaced seeds.

  # Helper to build a thin seed from current T
  def _thin_seed(current_T, max_take):
    n = len(current_T)
    if n == 0 or max_take <= 0:
      return []
    want = max(8, min(max_take, 40))
    step = max(1, n // want)
    return current_T[::step][:want]

  def build_micro_wave(current_T, budget, windows, connectors=True, alt_tag=0, reverse_blocks=False):
    if not current_T or budget <= 8:
      return []

    glo = min(l for l, _ in current_T)
    ghi = max(r for _, r in current_T)
    G = max(1, ghi - glo)

    # Thin, evenly spaced seed
    seed = _thin_seed(current_T, max(8, min(40, len(current_T) // 250)))
    if not seed:
      return []
    ulo = min(l for l, _ in seed)

    # Build translated micro-blocks at given windows; use staggered reversal to break symmetry
    blocks = []
    for j, (fa, fb) in enumerate(windows):
      win_lo = glo + int(round(fa * G))
      base = win_lo - ulo
      blk = [(l + base, r + base) for (l, r) in seed]
      if ((j + alt_tag) % 2) == 1:
        blk.reverse()
      blocks.append(blk)

    # Interleave micro-blocks (reverse global order when requested)
    micro = []
    maxlen = max(len(b) for b in blocks)
    order = list(range(len(blocks)))
    if reverse_blocks:
      order.reverse()
    for i in range(maxlen):
      for idx in order:
        b = blocks[idx]
        if i < len(b):
          micro.append(b[i])

    # Fractional-span connectors; moderate overlap to raise FF but control omega.
    if connectors:
      cons = [
        (glo + int(round(0.08 * G)), glo + int(round(0.30 * G))),  # left cap
        (glo + int(round(0.26 * G)), glo + int(round(0.56 * G))),  # mid cross
        (glo + int(round(0.44 * G)), glo + int(round(0.78 * G))),  # long cross
        (glo + int(round(0.60 * G)), glo + int(round(0.92 * G))),  # right cap
      ]
      for a, b in cons:
        if b > a:
          micro.append((a, b))

    return micro[:budget]

  # Phase 2A and 2B windows
  windows_A = [(0.12, 0.22), (0.35, 0.45), (0.58, 0.68), (0.80, 0.90)]
  windows_B = [(0.05, 0.15), (0.28, 0.38), (0.60, 0.70), (0.82, 0.92)]

  # Execute up to two micro-wave iterations driven by phase2_iters
  steps = min(max(0, int(phase2_iters)), 2)
  for iter_id in range(steps):
    room = CAP - len(T)
    if room <= 8:
      break
    # Alternate between A and B window families
    useA = (iter_id % 2 == 0)
    micro = build_micro_wave(
      T, room,
      windows_A if useA else windows_B,
      connectors=True,
      alt_tag=iter_id,
      reverse_blocks=(iter_id % 2 == 1),
    )
    if micro:
      T.extend(micro[:max(0, CAP - len(T))])

  # Secondary micro-wave: always attempt the complementary window family once
  room = CAP - len(T)
  if room > 8:
    # Complementary family with tiny edge pins and an extra long cross
    glo = min(l for l, _ in T)
    ghi = max(r for _, r in T)
    G = max(1, ghi - glo)
    micro_alt = build_micro_wave(T, room, windows_B, connectors=True, alt_tag=steps, reverse_blocks=True)
    # Extra long cross and tiny pins (capacity-guarded)
    extra = []
    extra.append((glo + int(round(0.18 * G)), glo + int(round(0.84 * G))))  # long cross
    eps = max(1, G // 512)
    for frac in (0.14, 0.86):
      a = glo + int(round(frac * G))
      b = a + eps
      if b > a:
        extra.append((a, b))
    combo = (micro_alt or []) + extra
    if combo:
      T.extend(combo[:max(0, CAP - len(T))])

  # Stage 3: Adaptive coverage-aware corridor injections
  if len(T) < CAP - 8:
    lo, hi, _ = _span_delta(T)
    spanG = max(1, hi - lo)

    # Coverage sampler to locate low-density regions
    def coverage_profile(seq, samples=64):
      if not seq:
        return [], []
      xs = []
      cnts = []
      # Sample strictly inside (open intervals): use midpoints of sub-buckets
      for i in range(samples):
        # Place sample point away from boundaries
        x = lo + 1 + int(((i + 0.5) / samples) * (spanG - 2))
        xs.append(x)
      # Count coverage: l < x < r (open intervals)
      for x in xs:
        c = 0
        for (l, r) in seq:
          if l < x < r:
            c += 1
        cnts.append(c)
      return xs, cnts

    xs, cnts = coverage_profile(T, samples=64)
    if xs:
      maxcov = max(cnts) if cnts else 0
      # Choose safe threshold: stay at least two below estimated max coverage
      thr_pin = max(0, min(maxcov - 2, 8))
      thr_conn = max(0, min(maxcov - 1, 9))

      safe_idxs = [i for i, c in enumerate(cnts) if c <= thr_pin]
      room = CAP - len(T)
      adaptive = []

      # Short pins at low-coverage locations (spread across span)
      if safe_idxs and room > 0:
        eps = max(1, spanG // 768)
        take = min(12, max(6, len(safe_idxs) // 6))
        step = max(1, len(safe_idxs) // take)
        for k in range(0, len(safe_idxs), step):
          i = safe_idxs[k]
          a = xs[i] - (eps // 2) + (k % 3)
          b = a + eps
          if b > a:
            adaptive.append((a, b))
          if len(adaptive) >= take:
            break

      # Mid-length corridor connectors along subsegments whose sampled coverage stays <= thr_conn
      if safe_idxs and (CAP - (len(T) + len(adaptive)) > 0):
        bridges = 0
        samples = len(xs)
        # Try a few candidates starting from leftmost safe indices
        for i in safe_idxs[:8]:
          j = min(samples - 1, i + samples // 8)  # ~12.5% of span length
          seg_max = max(cnts[i:j + 1]) if i < j else cnts[i]
          if seg_max <= thr_conn:
            L = xs[i]
            R = xs[j]
            if R > L:
              adaptive.append((L, R))
              bridges += 1
            if bridges >= 2:
              break

      # Capacity guard for adaptive phase
      if adaptive:
        room = CAP - len(T)
        if room > 0:
          T.extend(adaptive[:room])

  return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()