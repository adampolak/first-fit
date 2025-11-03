# EVOLVE-BLOCK-START

def construct_intervals(seed_count=1):
  """
  Construct a sequence of open intervals on the real line, in the order
  presented to FirstFit, to maximize FF colors divided by the clique number.

  Returns:
    intervals: list of (l, r) tuples (open intervals).
  """

  # Hard capacity guard
  CAP = 9800

  # Stable KT-spine start positions (Kierstead–Trotter pattern)
  spine_starts = (2, 6, 10, 14)

  # Seed: default single unit interval (best stability); multi-seed option retained
  if seed_count <= 1:
    T = [(0, 1)]
  else:
    step = 3
    T = [(i * step, i * step + 1) for i in range(max(1, int(seed_count)))]

  # -------------- Stage 1: KT spine, 6 rounds with mild interleaving --------------
  def apply_spine_round(current_T, ridx):
    lo = min(l for l, r in current_T)
    hi = max(r for l, r in current_T)
    delta = hi - lo
    if delta <= 0:
      delta = 1

    # Build four translated blocks
    blocks = []
    for start in spine_starts:
      base = start * delta - lo
      # Alternate inner orientation per (round, block) to enhance mixing deterministically
      block_T = current_T if ((ridx + start) % 2 == 0) else list(reversed(current_T))
      blocks.append([(l + base, r + base) for (l, r) in block_T])

    # Interleave on even rounds; sequential on odd rounds
    S = []
    if ridx % 2 == 0:
      maxlen = max(len(b) for b in blocks)
      for i in range(maxlen):
        for blk in blocks:
          if i < len(blk):
            S.append(blk[i])
    else:
      for blk in blocks:
        S.extend(blk)

    # Classic connectors (Figure 4 style) at delta scale
    s0, s1, s2, s3 = spine_starts
    connectors = [
      ((s0 - 1) * delta, (s1 - 1) * delta),  # left cap
      ((s2 + 2) * delta, (s3 + 2) * delta),  # right cap
      ((s0 + 2) * delta, (s2 - 1) * delta),  # cross 1
      ((s1 + 2) * delta, (s3 - 1) * delta),  # cross 2
    ]
    S.extend(connectors)
    return S

  # Run up to six rounds (size progression: 1 -> 8 -> 36 -> 148 -> 596 -> 2388 -> 9556)
  for ridx in range(6):
    nxt_size = 4 * len(T) + 4
    if nxt_size > CAP:
      break
    T = apply_spine_round(T, ridx)

  # If close to CAP, return backbone only
  if len(T) >= CAP - 16:
    return T

  # ---------- Utilities shared by micro-phases ----------
  def span_lo_hi(intervals):
    glo = min(l for l, r in intervals)
    ghi = max(r for l, r in intervals)
    return glo, ghi, max(1, ghi - glo)

  def thin_seed_even(current_T, max_seed):
    """Evenly spaced down-sample up to max_seed (deterministic)."""
    n = len(current_T)
    if n == 0 or max_seed <= 0:
      return []
    step = max(1, n // max_seed)
    return current_T[::step][:max_seed]

  def thin_seed_with_offset(current_T, max_seed, offset):
    """Evenly spaced sample with deterministic offset modulo the step."""
    n = len(current_T)
    if n == 0 or max_seed <= 0:
      return []
    step = max(1, n // max_seed)
    off = offset % step
    return current_T[off::step][:max_seed]

  def add_cross4_connectors(micro, delta2, glo, gate=4):
    """Add a small number of cross4 connectors: s0/s1 to s3 with safe offsets."""
    s0, s1, s2, s3 = spine_starts
    conns = [
      ((s0 + 1) * delta2 + glo, (s3 + 1) * delta2 + glo),
      ((s1 + 1) * delta2 + glo, (s3 + 2) * delta2 + glo),
      ((s0 + 2) * delta2 + glo, (s3 + 0) * delta2 + glo),
      ((s1 + 0) * delta2 + glo, (s3 + 3) * delta2 + glo),
    ]
    for a, b in conns[:max(0, gate)]:
      if b > a:
        micro.append((a, b))

  # ----------------- Stage 2a: Primary micro-phase (delta2 half-scale) -----------------
  def micro_round_primary(current_T, round_id, budget):
    if budget <= 0 or not current_T:
      return []

    glo, ghi, G = span_lo_hi(current_T)
    delta2 = max(1, G // 2)

    # Thin seed: bounded and deterministic
    # Slightly larger per-block target than previous to lift FF, still safe
    per_block_target = max(12, min(72, budget // 10))
    U = thin_seed_even(current_T, per_block_target)
    if not U:
      return []

    # Build four translated blocks with parity-based internal reversals
    blocks = []
    ulo = min(l for l, r in U)
    for s in spine_starts:
      base = s * delta2 - ulo
      block = [(l + base, r + base) for (l, r) in U]
      if ((s // 2) % 2) == (round_id % 2):
        block = list(reversed(block))
      blocks.append(block)

    # Interleave (forward on even round_id, reverse on odd)
    micro = []
    maxlen = max(len(b) for b in blocks)
    if round_id % 2 == 0:
      for i in range(maxlen):
        for blk in blocks:
          if i < len(blk):
            micro.append(blk[i])
    else:
      for i in range(maxlen):
        for blk in reversed(blocks):
          if i < len(blk):
            micro.append(blk[i])

    # Deterministic connectors at delta2 scale (classic + longer cross3)
    s0, s1, s2, s3 = spine_starts
    connectors = [
      ((s0 - 1) * delta2 + glo, (s1 - 1) * delta2 + glo),
      ((s2 + 2) * delta2 + glo, (s3 + 2) * delta2 + glo),
      ((s0 + 2) * delta2 + glo, (s2 - 1) * delta2 + glo),
      ((s1 + 2) * delta2 + glo, (s3 - 1) * delta2 + glo),
      ((s0 + 3) * delta2 + glo, (s3 + 3) * delta2 + glo),
    ]
    micro.extend(conn for conn in connectors if conn[1] > conn[0])

    # Add a tiny set of cross4 connectors to increase long-range coupling, gated
    add_cross4_connectors(micro, delta2, glo, gate=2)

    # Sparse micro caps: keep short and localized
    cap1 = (glo + (delta2 // 5), glo + int(1.6 * delta2))
    cap2 = (glo + int(0.95 * delta2), glo + int(2.4 * delta2))
    mid = glo + G // 2
    cap3 = (mid - max(1, delta2 // 10), mid + max(1, delta2 // 10))
    for cap in (cap1, cap2, cap3):
      if cap[1] > cap[0]:
        micro.append(cap)

    # Enforce budget
    if len(micro) > budget:
      micro = micro[:budget]
    return micro

  # ----------------- Stage 2b: Secondary window micro-phase (scaled embedding) -----------------
  # Window sets expressed in percent of global span; disjoint windows avoid clique spikes
  window_fracs_primary = [(7, 23), (31, 47), (55, 71), (79, 95)]
  window_fracs_secondary = [(5, 15), (28, 38), (60, 70), (82, 92)]

  def embed_scaled(U, a_pct, b_pct, glo, G):
    """Embed U linearly scaled to fit inside [glo + a%*G, glo + b%*G]."""
    width = max(1, (b_pct - a_pct) * G // 100)
    base = glo + (a_pct * G) // 100
    u_lo = min(l for l, r in U)
    u_hi = max(r for l, r in U)
    u_span = max(1, u_hi - u_lo)
    out = []
    for (l, r) in U:
      nl = base + ((l - u_lo) * width) // u_span
      nr = base + ((r - u_lo) * width) // u_span
      if nr <= nl:
        nr = nl + 1
      out.append((nl, nr))
    return out

  def window_micro_phase(current_T, round_id, budget, windows, seed_variation):
    """Build scaled copies of a thin seed inside designated windows; interleave across windows."""
    if budget <= 0 or not current_T:
      return []

    glo, ghi, G = span_lo_hi(current_T)

    # Deterministic seed derivation: Knuth constant mixer
    seed2 = (len(current_T) * 2654435761 + seed_variation) & 0xFFFFFFFF

    # Slightly larger target to ensure the windows are populated
    per_block_target = max(12, min(96, budget // 8))
    U = thin_seed_with_offset(current_T, per_block_target, offset=seed2)
    if not U:
      return []

    # Build one scaled block per window
    blocks = []
    for idx, (a, b) in enumerate(windows):
      block_U = list(reversed(U)) if ((round_id + idx) % 2 == 1) else U
      blk = embed_scaled(block_U, a, b, glo, G)
      blocks.append(blk)

    # Parity interleave across windows: forward if round_id even; reverse otherwise
    micro = []
    maxlen = max(len(b) for b in blocks)
    choose_blocks = blocks if (round_id % 2 == 0) else list(reversed(blocks))
    for i in range(maxlen):
      for blk in choose_blocks:
        if i < len(blk):
          micro.append(blk[i])

    # Short bridges between adjacent windows (localized to avoid omega blow-up)
    for j in range(len(windows) - 1):
      a1, b1 = windows[j]
      a2, b2 = windows[j + 1]
      left = glo + ((b1 - 2) * G) // 100
      right = glo + ((a2 + 2) * G) // 100
      if right > left:
        micro.append((left, right))

    # Add a very small number of delta2-scale cross4 connectors anchored near windows
    delta2 = max(1, G // 2)
    add_cross4_connectors(micro, delta2, glo, gate=2)

    # Enforce budget
    if len(micro) > budget:
      micro = micro[:budget]
    return micro

  # ----------------- Execute micro-phases under CAP -----------------
  remaining = CAP - len(T)
  # Primary micro-phase (delta2 half-scale)
  mr1_budget = max(0, remaining // 2)
  mr1 = micro_round_primary(T, round_id=0, budget=mr1_budget)
  if mr1:
    room = CAP - len(T)
    if len(mr1) > room:
      mr1 = mr1[:room]
    T.extend(mr1)

  # Secondary window micro-phase using two window sets back-to-back if room allows
  remaining = CAP - len(T)
  if remaining > 0:
    # First window set
    mr2_budget = max(0, remaining // 2)
    mr2 = window_micro_phase(T, round_id=1, budget=mr2_budget,
                             windows=window_fracs_primary, seed_variation=1)
    if mr2:
      room = CAP - len(T)
      if len(mr2) > room:
        mr2 = mr2[:room]
      T.extend(mr2)

  remaining = CAP - len(T)
  if remaining > 0:
    # Second window set (distinct) to increase late-stage diversity
    mr3_budget = remaining
    mr3 = window_micro_phase(T, round_id=2, budget=mr3_budget,
                             windows=window_fracs_secondary, seed_variation=2)
    if mr3:
      room = CAP - len(T)
      if len(mr3) > room:
        mr3 = mr3[:room]
      T.extend(mr3)

  return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()