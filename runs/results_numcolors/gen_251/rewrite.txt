# EVOLVE-BLOCK-START

def construct_intervals(seed_count=1):
  """
  Construct a sequence of intervals on the real line, in the order presented to FirstFit,
  engineered to drive up FirstFit's color usage while keeping the maximum clique (omega) modest.

  Interface compatibility:
    - Accepts seed_count (unused by evaluator, but supported).
    - Returns a list of (l, r) integer tuples with r > l (open intervals).
  """
  CAP = 9800
  CAP_MARGIN = 48  # headroom reserved for micro-phases

  # Rotating start templates (Kierstead–Trotter style)
  template_bank = [
    (2, 6, 10, 14),  # classic KT
    (1, 5, 9, 13),   # left-shifted
    (3, 7, 11, 15),  # right-shifted
    (4, 8, 12, 16),  # stretched-right
  ]

  # Deterministic utility: span stats
  def _span_delta(current_T):
    lo = min(l for l, r in current_T)
    hi = max(r for l, r in current_T)
    delta = hi - lo
    if delta <= 0:
      delta = 1
    return lo, hi, delta

  # Classic four KT connectors
  def _append_connectors(S, starts, delta):
    s0, s1, s2, s3 = starts
    pairs = [
      ((s0 - 1) * delta, (s1 - 1) * delta),  # left cap
      ((s2 + 2) * delta, (s3 + 2) * delta),  # right cap
      ((s0 + 2) * delta, (s2 - 1) * delta),  # cross 1
      ((s1 + 2) * delta, (s3 - 1) * delta),  # cross 2
    ]
    for a, b in pairs:
      if b > a:
        S.append((a, b))

  # Additional cross-scale connectors that couple both templates in a dual round
  def _append_union_connectors(S, startsA, startsB, delta):
    s_min = min(min(startsA), min(startsB))
    s_max = max(max(startsA), max(startsB))
    union_pairs = [
      ((s_min - 1) * delta, (s_max + 2) * delta),
      ((s_min + 2) * delta, (s_max - 1) * delta),
    ]
    for a, b in union_pairs:
      if b > a:
        S.append((a, b))
    # Cross-template mid-span couplers
    pairs2 = [
      ((startsA[1] + 1) * delta, (startsB[2] - 1) * delta),
      ((startsB[1] + 1) * delta, (startsA[2] - 1) * delta),
    ]
    for a, b in pairs2:
      if b > a:
        S.append((a, b))

  # Build four translated blocks for a given start template
  def _build_blocks(current_T, starts, lo, delta, parity_tag=0):
    blocks = []
    for idx, s in enumerate(starts):
      base = s * delta - lo
      # Alternate reversal by start index and tag to break symmetry
      if ((idx + parity_tag) % 2) == 1:
        src = list(reversed(current_T))
      else:
        src = current_T
      block = [(l + base, r + base) for (l, r) in src]
      blocks.append(block)
    return blocks

  # Interleave two templates' blocks: A0,B0,A1,B1,A2,B2,A3,B3 (or reversed on odd rounds)
  def _interleave_dual_blocks(blocksA, blocksB, reverse=False):
    order = []
    idxs = list(range(4))
    if reverse:
      idxs = list(reversed(idxs))
    for j in idxs:
      order.append(('A', j))
      order.append(('B', j))
    S = []
    m = max(len(blocksA[0]), len(blocksB[0])) if blocksA and blocksB else 0
    for i in range(m):
      for flag, j in order:
        blk = blocksA[j] if flag == 'A' else blocksB[j]
        if i < len(blk):
          S.append(blk[i])
    return S

  # One dual-template round (K=2) with deterministic interleaving and connectors
  def _dual_round(current_T, startsA, startsB, ridx):
    lo, hi, delta = _span_delta(current_T)
    # Parity tags diversify inner reversals; reverse interleaving on odd rounds
    blocksA = _build_blocks(current_T, startsA, lo, delta, parity_tag=(ridx % 2))
    blocksB = _build_blocks(current_T, startsB, lo, delta, parity_tag=((ridx + 1) % 2))
    S = _interleave_dual_blocks(blocksA, blocksB, reverse=(ridx % 2 == 1))
    # Append per-template and union connectors
    _append_connectors(S, startsA, delta)
    _append_connectors(S, startsB, delta)
    _append_union_connectors(S, startsA, startsB, delta)
    return S

  # Thin evenly-spaced seed from T
  def _thin_seed(current_T, max_seed):
    n = len(current_T)
    if n == 0 or max_seed <= 0:
      return []
    step = max(1, n // max_seed)
    return current_T[::step][:max_seed]

  # Four-window micro-round builder with deterministic connectors
  def _build_micro_round(current_T, window_fracs, budget, iter_id=0):
    if budget <= 8 or not current_T:
      return []
    glo = min(l for l, r in current_T)
    ghi = max(r for l, r in current_T)
    G = max(1, ghi - glo)

    # Thin seed for micro round (slightly larger on later iterations)
    seed_sz = max(12, min(40, len(current_T) // (280 if iter_id == 0 else 250)))
    U = _thin_seed(current_T, seed_sz)
    if not U:
      return []
    ulo = min(l for l, r in U)

    # Build translated micro-blocks aligned to windows, with small parity-based reversal
    blocks = []
    for k, (fa, fb) in enumerate(window_fracs):
      fa = max(0.05, min(0.90, fa))
      fb = max(0.10, min(0.95, fb))
      win_lo = glo + int(round(fa * G))
      base = win_lo - ulo
      block = [(l + base, r + base) for (l, r) in U]
      if ((k + iter_id) % 2) == 1:
        block = list(reversed(block))
      blocks.append(block)

    # Interleave the micro-blocks (reverse order on odd iterations)
    micro = []
    order = list(range(len(blocks)))
    if iter_id % 2 == 1:
      order.reverse()
    maxlen = max(len(b) for b in blocks)
    for i in range(maxlen):
      for idx in order:
        blk = blocks[idx]
        if i < len(blk):
          micro.append(blk[i])

    # Fractional-span connectors (analogous to KT caps)
    micro_connectors = [
      (glo + int(round(0.08 * G)), glo + int(round(0.30 * G))),
      (glo + int(round(0.60 * G)), glo + int(round(0.92 * G))),
      (glo + int(round(0.26 * G)), glo + int(round(0.56 * G))),
      (glo + int(round(0.44 * G)), glo + int(round(0.78 * G))),
      # One long-range connector to couple first and last windows
      (glo + int(round(0.12 * G)), glo + int(round(0.88 * G))),
    ]
    for a, b in micro_connectors:
      if b > a:
        micro.append((a, b))

    if len(micro) > budget:
      micro = micro[:budget]
    return micro

  # Insert a few long caps near the tail (late mixing with low clique impact)
  def _insert_tail_caps(T):
    lo, hi, delta = _span_delta(T)
    def cap_at(a_frac, b_frac):
      L = lo + max(1, int(round(a_frac * delta)))
      R = lo + max(1, int(round(b_frac * delta)))
      if R <= L:
        R = L + 1
      return (L, R)
    caps = [cap_at(0.10, 0.42), cap_at(0.36, 0.70), cap_at(0.62, 0.90)]
    out = list(T)
    for i, iv in enumerate(caps):
      if len(out) >= CAP:
        break
      pos = len(out) - (i * 3 + 1)
      if pos < 0:
        out.append(iv)
      else:
        out.insert(pos, iv)
    return out

  # Seed: single unit interval (supporting seed_count for compatibility)
  if seed_count <= 1:
    T = [(0, 1)]
  else:
    step = 3
    seeds = min(4, max(1, int(seed_count)))
    T = [(i * step, i * step + 1) for i in range(seeds)]

  # Stage 1: Dual-template KT spine for up to 4 rounds (K=2 density)
  # Pair templates deterministically: (0,1), (2,3), (1,2), (3,0)
  pairs = [(0, 1), (2, 3), (1, 2), (3, 0)]
  for ridx, (ia, ib) in enumerate(pairs):
    # Conservative capacity check (predictive): next size ~ 8*|T| + 12 connectors
    pred_next = 8 * len(T) + 16
    if pred_next > CAP - CAP_MARGIN:
      break
    startsA = template_bank[ia]
    startsB = template_bank[ib]
    T = _dual_round(T, startsA, startsB, ridx)
    if len(T) >= CAP - CAP_MARGIN:
      break

  # Early return if close to capacity
  if len(T) >= CAP - CAP_MARGIN:
    T = T[:CAP]
    # Normalize to non-negative integers and return
    min_l = min(l for l, r in T)
    if min_l < 0:
      T = [(l - min_l, r - min_l) for (l, r) in T]
    out = []
    for (l, r) in T:
      li = int(l)
      ri = int(r)
      if ri <= li:
        ri = li + 1
      out.append((li, ri))
    return out[:CAP]

  # Stage M1: First four-window micro-phase (gated)
  room = CAP - len(T)
  if room > CAP_MARGIN:
    windows1 = [(0.12, 0.22), (0.35, 0.45), (0.58, 0.68), (0.80, 0.90)]
    budget1 = min(room - CAP_MARGIN, 320)
    micro1 = _build_micro_round(T, windows1, budget=budget1, iter_id=0)
    if micro1:
      T.extend(micro1)

  # Stage M2: Second four-window micro-phase with distinct windows (gated)
  room = CAP - len(T)
  if room > CAP_MARGIN:
    windows2 = [(0.05, 0.15), (0.28, 0.38), (0.60, 0.70), (0.82, 0.92)]
    budget2 = min(room - CAP_MARGIN, 320)
    micro2 = _build_micro_round(T, windows2, budget=budget2, iter_id=1)
    if micro2:
      T.extend(micro2)

  # Tail caps
  if len(T) < CAP - 8:
    T = _insert_tail_caps(T)

  # Final capacity guard
  if len(T) > CAP:
    T = T[:CAP]

  # Normalize to non-negative integers and ensure r > l
  if not T:
    return []
  min_l = min(l for l, r in T)
  if min_l < 0:
    T = [(l - min_l, r - min_l) for (l, r) in T]
  intervals = []
  for (l, r) in T:
    li = int(l)
    ri = int(r)
    if ri <= li:
      ri = li + 1
    intervals.append((li, ri))
  return intervals[:CAP]

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()