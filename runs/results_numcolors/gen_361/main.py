# EVOLVE-BLOCK-START

def construct_intervals(seed_count=1):
  """
  Deterministic six-template spine plus a dual-pass micro-phase and tail connectors.
  Returns a list of (l, r) integer tuples (open intervals) in FirstFit presentation order.

  The design:
    - Six KT-like rounds using six distinct templates (one per round).
    - Deterministic interleaving/block-order schedules across rounds.
    - Small per-round pin connectors (short intervals) to couple colors safely.
    - Post-spine long-range connectors inserted near the tail.
    - Two CAP-aware micro passes with disjoint window families and connectors.

  Parameters:
    seed_count (int): supported for compatibility; uses 1 by default.

  Output:
    intervals: list[(int l, int r)], with r > l.
  """

  CAP = 9800

  # Six-template deterministic backbone (one per round)
  templates = [
    (2, 6, 10, 14),  # T0 classic KT
    (1, 5, 9, 13),   # T1 left-shifted
    (3, 7, 11, 15),  # T2 right-shifted
    (4, 8, 12, 16),  # T3 stretched-right
    (2, 5, 11, 14),  # T4 skewed symmetric
    (1, 7, 9, 15),   # T5 wide skew
  ]

  # Deterministic interleaving policy per round (True => interleave, False => sequential)
  interleave_schedule = [True, False, True, True, False, True]

  # Deterministic block order per round (indexes 0..3 into the four translated blocks)
  block_orders = [
    [0, 1, 2, 3],  # r0
    [2, 0, 3, 1],  # r1
    [1, 3, 0, 2],  # r2
    [3, 1, 2, 0],  # r3
    [0, 2, 1, 3],  # r4
    [2, 3, 1, 0],  # r5
  ]

  # Small jitter per round per block to break translation symmetry (kept modest to protect omega)
  jitter_sets = [
    (-0.020,  0.015, -0.010,  0.020),  # for r0
    ( 0.010, -0.015,  0.020, -0.010),  # for r1
    (-0.015,  0.010,  0.015, -0.005),  # for r2
    ( 0.005,  0.000, -0.005,  0.010),  # for r3
    ( 0.000,  0.005,  0.010, -0.015),  # for r4
    ( 0.015, -0.010,  0.005,  0.000),  # for r5
  ]

  # Seed with a single unit interval; support seed_count for compatibility
  if seed_count <= 1:
    T = [(0, 1)]
  else:
    seeds = min(4, max(1, int(seed_count)))
    step = 3
    T = [(i * step, i * step + 1) for i in range(seeds)]

  def span_delta(current_T):
    lo = min(l for l, r in current_T)
    hi = max(r for l, r in current_T)
    delta = hi - lo
    if delta <= 0:
      delta = 1
    return lo, hi, delta

  def append_classic_connectors(S, starts, delta):
    # KT-style four connectors
    s0, s1, s2, s3 = starts
    S.append(((s0 - 1) * delta, (s1 - 1) * delta))  # left cap
    S.append(((s2 + 2) * delta, (s3 + 2) * delta))  # right cap
    S.append(((s0 + 2) * delta, (s2 - 1) * delta))  # cross 1
    S.append(((s1 + 2) * delta, (s3 - 1) * delta))  # cross 2

  def add_round_pins(S, lo, delta, ridx):
    # Short pins to couple colors locally without raising omega much.
    # Enabled only on odd rounds to avoid stacking.
    if ridx % 2 == 1:
      eps = max(1, delta // 128)
      pin_positions = [0.18, 0.82]
      for pf in pin_positions:
        L = lo + int(round(pf * delta))
        R = L + eps
        if R > L:
          S.append((L, R))

  def apply_backbone_round(current_T, starts, ridx):
    lo, hi, delta = span_delta(current_T)
    jvec = jitter_sets[ridx]

    # Build four translated blocks with a small jitter perturbation
    blocks = []
    for idx, s in enumerate(starts):
      base = int(round((s + jvec[idx]) * delta)) - lo
      # reverse source ordering for odd blocks to reduce symmetry
      src = current_T[::-1] if (idx % 2 == 1) else current_T
      block = [(l + base, r + base) for (l, r) in src]
      blocks.append(block)

    # Compose S: either interleaved or sequential via a deterministic order
    S = []
    order = block_orders[ridx]
    if interleave_schedule[ridx]:
      maxlen = max(len(b) for b in blocks)
      for i in range(maxlen):
        for idx in order:
          blk = blocks[idx]
          if i < len(blk):
            S.append(blk[i])
    else:
      for idx in order:
        S.extend(blocks[idx])

    # Append classic connectors (scaled by delta)
    append_classic_connectors(S, starts, delta)
    # Add short pins for odd rounds (small, local coupling)
    add_round_pins(S, lo, delta, ridx)
    return S

  # Stage 1: Six-round deterministic spine
  for ridx in range(6):
    starts = templates[ridx]
    T = apply_backbone_round(T, starts, ridx)
    if len(T) >= CAP:
      T = T[:CAP]
      break

  if len(T) >= CAP - 8:
    # Near capacity; return strong baseline
    # Normalize and exit.
    intervals = []
    min_l = min(l for l, r in T) if T else 0
    if min_l < 0:
      T = [(l - min_l, r - min_l) for (l, r) in T]
    for (l, r) in T[:CAP]:
      li = int(l)
      ri = int(r)
      if ri <= li:
        ri = li + 1
      intervals.append((li, ri))
    return intervals

  # Utility: insert intervals near tail of sequence
  def insert_near_tail(seq, intervals):
    out = list(seq)
    for i, iv in enumerate(intervals):
      pos = len(out) - (2 * i + 1)
      if pos < 0:
        out.append(iv)
      else:
        out.insert(pos, iv)
    return out

  # Stage 1.5: deterministic long-range connectors added near the tail (post-spine only)
  lo, hi, G = span_delta(T)
  tail_connectors = [
    (lo + int(round(0.11 * G)), lo + int(round(0.51 * G))),
    (lo + int(round(0.24 * G)), lo + int(round(0.76 * G))),
    (lo + int(round(0.42 * G)), lo + int(round(0.86 * G))),
    (lo + int(round(0.18 * G)), lo + int(round(0.90 * G))),
  ]
  tail_connectors = [(a, b) for (a, b) in tail_connectors if b > a]
  room = CAP - len(T)
  if room > 0 and tail_connectors:
    T = insert_near_tail(T, tail_connectors[:room])

  if len(T) >= CAP - 8:
    intervals = []
    min_l = min(l for l, r in T) if T else 0
    if min_l < 0:
      T = [(l - min_l, r - min_l) for (l, r) in T]
    for (l, r) in T[:CAP]:
      li = int(l)
      ri = int(r)
      if ri <= li:
        ri = li + 1
      intervals.append((li, ri))
    return intervals

  # Thin evenly-spaced seed
  def thin_seed(seq, max_seed):
    n = len(seq)
    if n <= 0 or max_seed <= 0:
      return []
    step = max(1, n // max_seed)
    return seq[::step][:max_seed]

  # Build a micro pass: translate a thin seed to multiple windows and interleave
  def build_micro_pass(current_T, budget, windows, seed_cap, add_long_caps=False, add_pins=True):
    if not current_T or budget <= 8:
      return []

    glo = min(l for l, r in current_T)
    ghi = max(r for l, r in current_T)
    Gspan = max(1, ghi - glo)

    # Thin seed
    seed_sz = max(8, min(seed_cap, max(8, len(current_T) // 280)))
    U = thin_seed(current_T, seed_sz)
    if not U:
      return []

    ulo = min(l for l, r in U)

    # Translate into windows
    blocks = []
    for (fa, fb) in windows:
      win_lo = glo + int(round(fa * Gspan))
      base = win_lo - ulo
      blk = [(l + base, r + base) for (l, r) in U]
      # reverse internally for half of the windows to break symmetry
      if int(round(fa * 100)) // 5 % 2 == 1:
        blk = list(reversed(blk))
      blocks.append(blk)

    # Interleave blocks deterministically
    micro = []
    maxlen = max(len(b) for b in blocks)
    # alternate order pattern
    order = list(range(len(blocks)))
    if len(blocks) >= 4:
      order = [0, 2, 1, 3] + list(range(4, len(blocks)))
    for i in range(maxlen):
      for idx in order:
        blk = blocks[idx]
        if i < len(blk):
          micro.append(blk[i])

    # Fractional-span connectors (lightweight) across windows
    connectors = [
      (glo + int(round(0.08 * Gspan)), glo + int(round(0.30 * Gspan))),
      (glo + int(round(0.62 * Gspan)), glo + int(round(0.94 * Gspan))),
      (glo + int(round(0.26 * Gspan)), glo + int(round(0.56 * Gspan))),
      (glo + int(round(0.44 * Gspan)), glo + int(round(0.78 * Gspan))),
    ]
    if add_long_caps:
      # One extra long connector gated to avoid too much span coverage
      connectors.append((glo + int(round(0.20 * Gspan)), glo + int(round(0.84 * Gspan))))
    for (a, b) in connectors:
      if b > a:
        micro.append((a, b))

    # Short pins to lightly couple without omega blowup
    if add_pins:
      eps = max(1, Gspan // 512)
      for pf in (0.16, 0.39, 0.61, 0.85):
        L = glo + int(round(pf * Gspan))
        R = L + eps
        if R > L:
          micro.append((L, R))

    if len(micro) > budget:
      micro = micro[:budget]
    return micro

  # Stage 2: two micro passes, CAP-aware, disjoint windows, deterministic
  room = CAP - len(T)

  if room > 8:
    # Pass A windows (inner)
    windows_A = [(0.10, 0.18), (0.32, 0.40), (0.54, 0.62), (0.76, 0.84)]
    budget_A = max(0, room // 2)
    micro_A = build_micro_pass(T, budget_A, windows_A, seed_cap=32, add_long_caps=True, add_pins=True)
    if micro_A:
      avail = CAP - len(T)
      if len(micro_A) > avail:
        micro_A = micro_A[:avail]
      T.extend(micro_A)

  room = CAP - len(T)
  if room > 8:
    # Pass B windows (offset, outer)
    windows_B = [(0.18, 0.26), (0.40, 0.48), (0.62, 0.70), (0.84, 0.92)]
    budget_B = max(0, room)
    micro_B = build_micro_pass(T, budget_B, windows_B, seed_cap=40, add_long_caps=False, add_pins=True)
    if micro_B:
      avail = CAP - len(T)
      if len(micro_B) > avail:
        micro_B = micro_B[:avail]
      T.extend(micro_B)

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
  for (l, r) in T[:CAP]:
    li = int(l)
    ri = int(r)
    if ri <= li:
      ri = li + 1
    intervals.append((li, ri))

  return intervals

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()