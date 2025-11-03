# EVOLVE-BLOCK-START

def construct_intervals(enable_alt_microphase=True):
  """
  Planner-guided six-template backbone with two CAP-aware micro-phases
  and per-round density boosters. Returns a list of (l, r) integer pairs.
  """

  # Hard capacity guard (< 10000)
  CAP = 9800
  BASE_SEED = 137  # deterministic driver for template selection and interleaving

  # Six deterministic start-pattern templates (4, 5, and 6 blocks)
  templates = [
    (2, 6, 10, 14),       # T0 classic KT
    (1, 5, 9, 13),        # T1 left-shifted
    (3, 7, 11, 15),       # T2 right-shifted
    (4, 8, 12, 16),       # T3 stretched-right
    (2, 5, 8, 11, 14),    # T4 5-block variant
    (2, 4, 7, 10, 13, 16) # T5 6-block variant
  ]

  # Seed
  T = [(0, 1)]

  # ------------------ Helpers ------------------

  def _span(Tcur):
    lo = min(l for l, _ in Tcur)
    hi = max(r for _, r in Tcur)
    d = hi - lo
    return lo, hi, (1 if d <= 0 else d)

  def _interleave_blocks(blocks, reverse=False):
    S = []
    maxlen = max((len(b) for b in blocks), default=0)
    order = list(range(len(blocks)))
    if reverse:
      order.reverse()
    for i in range(maxlen):
      for idx in order:
        blk = blocks[idx]
        if i < len(blk):
          S.append(blk[i])
    return S

  def _concat_blocks(blocks, reverse=False):
    if reverse:
      blocks = list(reversed(blocks))
    S = []
    for blk in blocks:
      S.extend(blk)
    return S

  def _general_connectors(starts, delta, limit=6):
    # Left / right caps plus sparse cross-links
    k = len(starts)
    cons = []
    if k >= 2:
      cons.append(((starts[0] - 1) * delta, (starts[1] - 1) * delta))          # left cap
      cons.append(((starts[-2] + 2) * delta, (starts[-1] + 2) * delta))        # right cap
    # Cross links: (s_i+2, s_{i+2}-1), sparsified to avoid omega spikes
    for i in range(max(0, k - 2)):
      a = (starts[i] + 2) * delta
      b = (starts[i + 2] - 1) * delta
      if b > a:
        cons.append((a, b))
      if len(cons) >= limit:
        break
    return cons[:limit]

  def _apply_template_round(Tcur, starts, interleave=True, reverse=False, connector_cap=6):
    lo, hi, delta = _span(Tcur)
    # Build translated blocks
    blocks = []
    for s in starts:
      base = s * delta - lo
      block = [(l + base, r + base) for (l, r) in Tcur]
      blocks.append(block)
    # Assemble
    if interleave:
      S = _interleave_blocks(blocks, reverse=reverse)
    else:
      S = _concat_blocks(blocks, reverse=reverse)
    # General connectors (sparse and bounded)
    S.extend(_general_connectors(starts, delta, limit=connector_cap))
    return S

  def _predict_size(sz, k, conn_cnt=6):
    # conservative predictor for next round size
    return k * sz + conn_cnt

  def _tail_caps(Tcur, caps, cap_limit):
    # Insert caps near the tail in a stable pattern
    out = list(Tcur)
    room = cap_limit - len(out)
    if room <= 0:
      return out
    for i, iv in enumerate(caps[:room]):
      pos = len(out) - (2 * i + 1)
      if pos <= 0:
        out.append(iv)
      else:
        out.insert(pos, iv)
    return out

  def _density_boost(Tcur, frac_windows, per_round_quota=24):
    # Per-round thin boost: compress a tiny seed into inner windows (very short intervals)
    if per_round_quota <= 0 or not Tcur:
      return []
    lo, hi, d = _span(Tcur)
    G = d
    # Thin seed (very small)
    seed_sz = max(8, min(16, len(Tcur) // 600))
    stride = max(1, len(Tcur) // max(1, seed_sz))
    U = [Tcur[i] for i in range(0, len(Tcur), stride)][:seed_sz]
    if not U:
      return []
    ulo = min(l for l, _ in U)
    eps = max(1, G // 1024)
    micro = []
    for (fa, fb) in frac_windows:
      win_lo = lo + int(fa * G)
      base = win_lo - ulo
      for idx, (l, r) in enumerate(U):
        mid = (l + r) // 2
        L = mid + base + (idx % 3)  # mild stagger
        R = L + eps
        if R > L:
          micro.append((L, R))
        if len(micro) >= per_round_quota:
          break
      if len(micro) >= per_round_quota:
        break
    return micro

  def _long_range_connectors(lo, hi):
    G = max(1, hi - lo)
    # Deterministic long connectors spanning across the backbone, carefully spaced
    L = [
      (lo + int(0.08 * G), lo + int(0.60 * G)),
      (lo + int(0.25 * G), lo + int(0.75 * G)),
      (lo + int(0.44 * G), lo + int(0.78 * G)),
      (lo + int(0.60 * G), lo + int(0.92 * G)),
      (lo + int(0.18 * G), lo + int(0.84 * G)),
      (lo + int(0.12 * G), lo + int(0.33 * G)),
    ]
    # ensure r>l strictly
    L2 = []
    for a, b in L:
      if b <= a:
        b = a + 1
      L2.append((a, b))
    return L2

  # ------------------ Stage 1: Planner-guided six-template backbone ------------------

  # Deterministic round budget: aim for 6 rounds but back off if CAP predicts overflow
  max_rounds = 6
  for ridx in range(max_rounds):
    # Choose template deterministically from BASE_SEED and round
    sel = (BASE_SEED + 97 * ridx) % len(templates)
    starts = templates[sel]
    k = len(starts)

    # Interleaving policy: enable on even rounds; reverse order on rounds with (ridx+BASE_SEED) odd
    do_inter = (ridx % 2 == 0)
    rev = ((ridx + BASE_SEED) % 2 == 1)

    # Predict size; if overflow, fall back to classic 4-block template
    next_size = _predict_size(len(T), k, conn_cnt=6)
    if next_size > CAP - 300:  # reserve micro-phase budget
      starts = templates[0]  # fallback to (2,6,10,14)
      k = len(starts)
      next_size = _predict_size(len(T), k, conn_cnt=6)
      if next_size > CAP - 300:
        # If still too big, stop backbone
        break

    # Apply round
    T = _apply_template_round(T, starts, interleave=do_inter, reverse=rev, connector_cap=6)

    # Per-round density boost (very thin to avoid omega spikes), deterministic windows
    if len(T) < CAP - 64:
      frac_windows = [(0.31, 0.33), (0.49, 0.51), (0.67, 0.69)]
      boost = _density_boost(T, frac_windows, per_round_quota=24)
      room = CAP - len(T)
      if boost and room > 0:
        if len(boost) > room:
          boost = boost[:room]
        T.extend(boost)

  # If nearly at capacity, return early
  if len(T) >= CAP - 12:
    return T[:CAP]

  # ------------------ Deterministic long-range connectors across the finished spine ------------------
  lo, hi, _ = _span(T)
  lr = _long_range_connectors(lo, hi)
  room = CAP - len(T)
  if room > 0 and lr:
    if len(lr) > room:
      lr = lr[:room]
    T.extend(lr)

  if len(T) >= CAP - 12:
    T = T[:CAP]
    return T

  # ------------------ Two-phase CAP-aware micro-phases ------------------

  def _thin_seed(Tcur, max_seed):
    n = len(Tcur)
    if n <= 0 or max_seed <= 0:
      return []
    step = max(1, n // max_seed)
    return Tcur[::step][:max_seed]

  def _build_micro_phase(Tcur, budget, config_id=0):
    # Comb/tower hybrid micro-phase with connectors; CAP-aware output size
    if budget <= 8 or not Tcur:
      return []

    glo, ghi, G = _span(Tcur)

    # Seed size by config
    base_seed = max(12, min(48, len(Tcur) // 220))
    if config_id == 1:
      base_seed = max(16, min(64, len(Tcur) // 180))
    U = _thin_seed(Tcur, base_seed)
    if not U:
      return []

    ulo = min(l for l, _ in U)

    # Window sets
    if config_id == 0:
      windows = [(0.10, 0.18), (0.32, 0.40), (0.54, 0.62), (0.76, 0.84)]
      order = [0, 2, 1, 3]
    else:
      windows = [(0.06, 0.14), (0.28, 0.36), (0.60, 0.68), (0.82, 0.90)]
      order = [3, 1, 0, 2]

    eps = max(1, G // (768 if config_id == 0 else 640))
    levels = 2 if config_id == 0 else 3

    # Build tower/comb blocks
    blocks = []
    for t_idx, (fa, fb) in enumerate(windows):
      win_lo = glo + int(round(fa * G))
      base = win_lo - ulo
      block = []
      for idx, (l, r) in enumerate(U):
        mid = (l + r) // 2
        for lv in range(levels):
          jitter = (idx + lv + 2 * t_idx) % 4
          L = mid + base - (eps // 2) + jitter + lv
          R = L + eps
          if R > L:
            block.append((L, R))
      blocks.append(block)

    # Interleave blocks with a deterministic order
    micro = []
    maxlen = max((len(b) for b in blocks), default=0)
    for i in range(maxlen):
      for idx in order:
        blk = blocks[idx]
        if i < len(blk):
          micro.append(blk[i])

    # Add fractional connectors to tie windows (avoid too many)
    cons = [
      (glo + int(round(0.09 * G)), glo + int(round(0.30 * G))),
      (glo + int(round(0.42 * G)), glo + int(round(0.74 * G))),
      (glo + int(round(0.64 * G)), glo + int(round(0.91 * G))),
    ]
    for a, b in cons:
      if b > a:
        micro.append((a, b))

    # Trim to budget
    if len(micro) > budget:
      micro = micro[:budget]
    return micro

  # Budgeting: split the remaining CAP into two micro phases (60% / 40%)
  remaining = CAP - len(T)
  if remaining > 16:
    bud1 = max(0, min(remaining - 8, int(0.60 * remaining)))
    bud2 = max(0, remaining - bud1)

    # Phase 1
    micro1 = _build_micro_phase(T, bud1, config_id=0)
    if micro1:
      avail = CAP - len(T)
      if len(micro1) > avail:
        micro1 = micro1[:avail]
      T.extend(micro1)

    # Phase 2 (alternate windows)
    if enable_alt_microphase:
      avail = CAP - len(T)
      if avail > 8:
        micro2 = _build_micro_phase(T, avail, config_id=1)
        if micro2:
          if len(micro2) > avail:
            micro2 = micro2[:avail]
          T.extend(micro2)

  # Final normalization and CAP trim
  if len(T) > CAP:
    T = T[:CAP]

  # Ensure integer endpoints and r > l
  if T:
    min_l = min(l for l, _ in T)
    if min_l < 0:
      T = [(l - min_l, r - min_l) for (l, r) in T]
  out = []
  for (l, r) in T:
    li = int(l)
    ri = int(r)
    if ri <= li:
      ri = li + 1
    out.append((li, ri))
  if len(out) > CAP:
    out = out[:CAP]
  return out

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()