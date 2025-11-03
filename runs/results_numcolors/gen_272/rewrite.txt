# EVOLVE-BLOCK-START

def construct_intervals(seed_count=1):
  """
  Deterministic KT-style spine with parity interleaving, dual micro-phases,
  and CAP-aware long-range connectors to raise FirstFit while keeping omega ~<= 10.

  Inputs:
    seed_count (int): retained for interface compatibility; nonnegative.
  Returns:
    List[(l, r)] of open intervals with integer endpoints and r > l.
  """
  # Hard capacity to satisfy evaluator constraint
  CAP = 9800

  # Deterministic base seed (no RNG used; just a constant for future-proofing)
  BASE_SEED = 1337

  # Start templates (Kierstead–Trotter inspired)
  template_bank_primary = [
    (2, 6, 10, 14),  # T1: classic KT
    (1, 5, 9, 13),   # T2: left-shifted
    (3, 7, 11, 15),  # T3: right-shifted
    (4, 8, 12, 16),  # T4: stretched-right
  ]
  # Auxiliary templates used selectively to add diversity without densifying too much
  template_bank_aux = [
    (2, 4, 8, 12),   # inner-left biased
    (3, 5, 9, 13),   # gentle left pack
    (1, 7, 11, 15),  # wide skew
  ]

  # Seed
  if seed_count <= 1:
    T = [(0, 1)]
  else:
    # Multiple disjoint seeds; keep small to avoid early omega jump
    seeds = min(4, max(1, int(seed_count)))
    step = 3
    T = [(i * step, i * step + 1) for i in range(seeds)]

  # Helpers
  def _span_delta(current_T):
    lo = min(l for l, r in current_T)
    hi = max(r for l, r in current_T)
    delta = hi - lo
    if delta <= 0:
      delta = 1
    return lo, hi, delta

  def _append_kt_connectors(S, starts, delta):
    s0, s1, s2, s3 = starts
    S.append(((s0 - 1) * delta, (s1 - 1) * delta))  # left cap
    S.append(((s2 + 2) * delta, (s3 + 2) * delta))  # right cap
    S.append(((s0 + 2) * delta, (s2 - 1) * delta))  # cross 1
    S.append(((s1 + 2) * delta, (s3 - 1) * delta))  # cross 2

  def _append_long_range_connectors(S, starts, delta, budget=2):
    # Deterministic long-range connectors spanning 4–6 delta to enhance cross-scale coupling
    # Anchored to the extremal starts. Keep small (<= 2) to avoid omega inflation.
    s0, s1, s2, s3 = starts
    cands = [
      ((s0 + 4) * delta, (s3 + 4) * delta),
      ((s0 + 5) * delta, (s3 + 6) * delta),
      ((s1 + 4) * delta, (s2 + 5) * delta),
    ]
    added = 0
    for a, b in cands:
      if b > a:
        S.append((a, b))
        added += 1
        if added >= max(0, int(budget)):
          break

  def _apply_round(current_T, starts, do_interleave=False, reverse_order=False, add_lr=False):
    lo, hi, delta = _span_delta(current_T)

    # Build four translated blocks
    blocks = []
    for s in starts:
      base = s * delta - lo
      block = [(l + base, r + base) for (l, r) in current_T]
      blocks.append(block)

    # Assemble S
    S = []
    if do_interleave:
      order = list(range(4))
      if reverse_order:
        order.reverse()
      maxlen = max(len(b) for b in blocks)
      for i in range(maxlen):
        for idx in order:
          blk = blocks[idx]
          if i < len(blk):
            S.append(blk[i])
    else:
      seq_blocks = list(reversed(blocks)) if reverse_order else blocks
      for blk in seq_blocks:
        S.extend(blk)

    # Append classic KT connectors
    _append_kt_connectors(S, starts, delta)

    # Optional long-range connectors (kept conservative)
    if add_lr:
      _append_long_range_connectors(S, starts, delta, budget=2)

    return S

  # Predictive size accounting: KT growth per round size -> 4*size + 4
  def _can_advance_round(sz):
    return (4 * sz + 4) <= CAP

  # Stage 1: spine assembly with deterministic parity interleaving
  round_idx = 0
  while _can_advance_round(len(T)) and round_idx < 7:
    # Primary templates rotate every round; every third round uses an auxiliary template deterministically
    if round_idx % 3 == 2:
      starts = template_bank_aux[(BASE_SEED + round_idx) % len(template_bank_aux)]
    else:
      starts = template_bank_primary[(BASE_SEED + round_idx) % len(template_bank_primary)]

    # Parity policy: interleave on even rounds; reverse order on odd rounds
    do_interleave = (round_idx % 2 == 0)
    reverse_order = (round_idx % 2 == 1)

    # Add long-range connectors only on rounds 2 and 5 (0-based), conservative gating
    add_lr = (round_idx in (2, 5))

    T = _apply_round(T, starts, do_interleave=do_interleave, reverse_order=reverse_order, add_lr=add_lr)
    if len(T) >= CAP:
      T = T[:CAP]
      break
    round_idx += 1

  # Early return if near capacity
  if len(T) >= CAP - 8:
    # Normalize to integer open intervals
    return _normalize_intervals(T, CAP)

  # Micro-phase A: fractional-window thin round (primary windows)
  def _thin_seed(current_T, max_seed):
    n = len(current_T)
    if n <= 0 or max_seed <= 0:
      return []
    step = max(1, n // max_seed)
    return current_T[::step][:max_seed]

  def _build_micro(current_T, budget, windows, seed_cap, alt_tag=0, add_lr=False):
    if not current_T or budget <= 8:
      return []

    glo = min(l for l, r in current_T)
    ghi = max(r for l, r in current_T)
    G = max(1, ghi - glo)

    # Thin seed
    seed_sz = max(8, min(seed_cap, max(8, len(current_T) // 280)))
    U = _thin_seed(current_T, seed_sz)
    if not U:
      return []
    ulo = min(l for l, r in U)

    # Blocks translated into windows
    blocks = []
    for (fa, fb) in windows:
      fa = max(0.05, min(0.90, fa))
      fb = max(0.10, min(0.95, fb))
      win_lo = glo + int(round(fa * G))
      base = win_lo - ulo
      block = [(l + base, r + base) for (l, r) in U]
      # Parity-breaking reversal
      tag = alt_tag + int(round(fa * 100)) // 5
      if tag % 2 == 1:
        block = list(reversed(block))
      blocks.append(block)

    # Interleave micro-blocks; reverse block order for odd alt_tag
    micro = []
    order = list(range(len(blocks)))
    if (alt_tag % 2) == 1:
      order.reverse()
    maxlen = max(len(b) for b in blocks)
    for i in range(maxlen):
      for idx in order:
        blk = blocks[idx]
        if i < len(blk):
          micro.append(blk[i])

    # Add fractional connectors (KT analog)
    micro_connectors = [
      (glo + int(round(0.08 * G)), glo + int(round(0.30 * G))),
      (glo + int(round(0.60 * G)), glo + int(round(0.92 * G))),
      (glo + int(round(0.26 * G)), glo + int(round(0.56 * G))),
      (glo + int(round(0.44 * G)), glo + int(round(0.78 * G))),
    ]
    for a, b in micro_connectors:
      if b > a:
        micro.append((a, b))

    # Optional long-range connector at the micro scale
    if add_lr:
      lr = (glo + int(round(0.18 * G)), glo + int(round(0.84 * G)))
      if lr[1] > lr[0]:
        micro.append(lr)

    # Trim by budget
    if len(micro) > budget:
      micro = micro[:budget]
    return micro

  # Windows for two deterministic micro-phases
  window_fracs_A = [(0.12, 0.22), (0.35, 0.45), (0.58, 0.68), (0.80, 0.90)]
  window_fracs_B = [(0.06, 0.14), (0.28, 0.38), (0.54, 0.64), (0.78, 0.88)]

  # Execute Micro A then Micro B, each once, capacity-guarded
  room = CAP - len(T)
  if room > 8:
    microA = _build_micro(T, room, window_fracs_A, seed_cap=36, alt_tag=(BASE_SEED % 2), add_lr=False)
    if microA:
      if len(microA) > room:
        microA = microA[:room]
      T.extend(microA)

  room = CAP - len(T)
  if room > 8:
    microB = _build_micro(T, room, window_fracs_B, seed_cap=40, alt_tag=((BASE_SEED + 1) % 4), add_lr=True)
    if microB:
      if len(microB) > room:
        microB = microB[:room]
      T.extend(microB)

  # Post-micro long-range connectors across the global span (very few, CAP-aware)
  if len(T) < CAP - 4:
    lo, hi, delta = _span_delta(T)
    # Deterministic cross-scale caps; 2 pieces, conservative
    post_caps = [
      (lo + max(1, int(round(0.15 * delta))), lo + max(1, int(round(0.60 * delta)))),
      (lo + max(1, int(round(0.40 * delta))), lo + max(1, int(round(0.85 * delta)))),
    ]
    add = []
    for a, b in post_caps:
      if b > a:
        add.append((a, b))
    room = CAP - len(T)
    add = add[:room]
    T.extend(add)

  # Final capacity guard and normalization
  return _normalize_intervals(T, CAP)


# Utility: normalize to integer open intervals and cap length
def _normalize_intervals(T, CAP):
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
  if len(intervals) > CAP:
    intervals = intervals[:CAP]
  return intervals

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()