# EVOLVE-BLOCK-START

def construct_intervals(enable_alt_microphase=True):
  """
  Deterministic seeded KT-style construction with CAP-aware densification and
  two micro-phases to increase FirstFit colors while controlling omega.

  Args:
    enable_alt_microphase (bool): whether to enable the alternate micro-phase.

  Returns:
    list of (l, r) tuples representing open intervals presented in order.
  """

  # -------------------------
  # Global controls and seed
  # -------------------------
  CAP = 9800                 # hard cap; keep < 10000 for evaluation headroom
  BASE_SEED = 91138233       # deterministic, centralized seed for all decisions
  MAX_ROUNDS = 6             # backbone rounds; KT-style growth under CAP

  # Start-template bank (rotated deterministically per round)
  TEMPLATE_BANK = [
    (2, 6, 10, 14),
    (1, 5, 9, 13),
    (3, 7, 11, 15),
    (4, 8, 12, 16),
  ]

  # -------------------------
  # Utility functions
  # -------------------------
  MASK64 = (1 << 64) - 1
  def _mix64(x):
    z = (x + 0x9E3779B97F4A7C15) & MASK64
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9 & MASK64
    z = (z ^ (z >> 27)) * 0x94D049BB133111EB & MASK64
    z ^= (z >> 31)
    return z & MASK64

  def _per_round_seed(base, ridx):
    return _mix64(base ^ _mix64(ridx + 0xD1B54A32D192ED03))

  def _span(intervals):
    lo = min(l for l, r in intervals)
    hi = max(r for l, r in intervals)
    delta = hi - lo
    if delta <= 0:
      delta = 1
    return lo, hi, delta

  def _translate_block(src, base_shift):
    # Translate src by base_shift; endpoints can be fractional; we normalize later.
    return [(l + base_shift, r + base_shift) for (l, r) in src]

  def _interleave_blocks(blocks, reverse=False, rotate=0):
    if not blocks:
      return []
    k = len(blocks)
    order = list(range(k))
    if reverse:
      order.reverse()
    if k > 0 and rotate % k != 0:
      rot = rotate % k
      order = order[rot:] + order[:rot]
    out = []
    maxlen = max(len(b) for b in blocks)
    for i in range(maxlen):
      for idx in order:
        blk = blocks[idx]
        if i < len(blk):
          out.append(blk[i])
    return out

  def _add_connectors(starts, delta, origin=0.0):
    s0, s1, s2, s3 = starts
    return [
      (origin + (s0 - 1) * delta, origin + (s1 - 1) * delta),  # left cap
      (origin + (s2 + 2) * delta, origin + (s3 + 2) * delta),  # right cap
      (origin + (s0 + 2) * delta, origin + (s2 - 1) * delta),  # cross 1
      (origin + (s1 + 2) * delta, origin + (s3 - 1) * delta),  # cross 2
    ]

  def _thin_seed(intervals, max_seed):
    n = len(intervals)
    if n == 0 or max_seed <= 0:
      return []
    step = max(1, n // max_seed)
    return intervals[::step][:max_seed]

  def _normalize(intervals):
    if not intervals:
      return []
    # Shift to non-negative if needed
    min_l = min(l for l, _ in intervals)
    if min_l < 0:
      shift = -min_l
      intervals = [(l + shift, r + shift) for (l, r) in intervals]
    # Cast to ints and ensure r > l
    norm = []
    for (l, r) in intervals:
      li = int(l)
      ri = int(r)
      if ri <= li:
        ri = li + 1
      norm.append((li, ri))
    return norm

  # Long-range deterministic connectors (4–6 delta spans) with conservative count
  def _long_range_connectors(lo, hi):
    G = max(1, hi - lo)
    return [
      (lo + int(0.06 * G), lo + int(0.40 * G)),
      (lo + int(0.60 * G), lo + int(0.94 * G)),
    ]

  # -------------------------
  # Backbone (KT-style) with parity interleaving and CAP-aware densification
  # -------------------------
  T = [(0, 1)]  # seed

  for ridx in range(MAX_ROUNDS):
    # Predict rough size after this round to keep within CAP: sz -> 4*sz + 4 (+ small densification)
    projected = 4 * len(T) + 4
    if projected > CAP:
      break

    lo, hi, delta = _span(T)
    starts = TEMPLATE_BANK[ridx % len(TEMPLATE_BANK)]

    # Build 4 translated blocks
    blocks = []
    for s in starts:
      base = s * delta - lo
      blocks.append(_translate_block(T, base))

    # Parity schedule: even rounds interleave forward; odd rounds interleave reversed
    do_reverse = (ridx % 2 == 1)
    S = _interleave_blocks(blocks, reverse=do_reverse, rotate=ridx)

    # Classic 4 connectors
    S.extend(_add_connectors(starts, delta, origin=0.0))

    # CAP-aware densification: add thin ghost blocks at half-offset positions to raise FF pressure
    # without significantly inflating omega (bounded thin seed, staggered placement).
    room = CAP - len(S)
    if room > 0:
      # small thin seed (size scales sublinearly with |T|)
      seed_sz = max(6, min(24, len(T) // 350))
      U = _thin_seed(T, seed_sz)
      if U:
        ulo = min(l for l, r in U)
        # Half-offset windows anchored between existing starts (use two ghosts)
        ghost_starts = [starts[0] + 0.5, starts[2] + 0.5]
        ghost_blocks = []
        for gs in ghost_starts:
          gbase = gs * delta - ulo
          ghost_blocks.append(_translate_block(U, gbase))
        ghosts = _interleave_blocks(ghost_blocks, reverse=False, rotate=ridx + 1)

        # Trim ghosts if low capacity
        if len(ghosts) > max(0, room // 2):
          ghosts = ghosts[:max(0, room // 2)]
        S.extend(ghosts)

    T = S

  # Early return if at capacity
  if len(T) >= CAP:
    return _normalize(T)

  # -------------------------
  # Micro-phase 1: fractional-window densification (delta-scale)
  # -------------------------
  def _micro_phase(Tin, budget, round_hint):
    if budget <= 0 or not Tin:
      return []
    glo = min(l for l, r in Tin)
    ghi = max(r for l, r in Tin)
    G = max(1, ghi - glo)

    # Thin seed; bounded to keep clique small
    seed_sz = max(8, min(32, len(Tin) // 300))
    U = _thin_seed(Tin, seed_sz)
    if not U:
      return []
    ulo = min(l for l, r in U)

    window_fracs = [(0.12, 0.22), (0.35, 0.45), (0.58, 0.68), (0.80, 0.90)]
    blocks = []
    for (fa, fb) in window_fracs:
      win_lo = glo + int(round(fa * G))
      base = win_lo - ulo
      blocks.append(_translate_block(U, base))

    # Interleave forward on even hint, reverse on odd hint
    micro = _interleave_blocks(blocks, reverse=(round_hint % 2 == 1), rotate=round_hint)

    # Micro-scale connectors within [glo, ghi]
    micro_connectors = [
      (glo + int(0.08 * G), glo + int(0.30 * G)),
      (glo + int(0.60 * G), glo + int(0.92 * G)),
      (glo + int(0.26 * G), glo + int(0.56 * G)),
      (glo + int(0.44 * G), glo + int(0.78 * G)),
    ]
    micro.extend(micro_connectors)

    # Budget trim
    if len(micro) > budget:
      micro = micro[:budget]
    return micro

  room = CAP - len(T)
  if room > 8:
    micro1 = _micro_phase(T, room, round_hint=0)
    if micro1:
      if len(micro1) > room:
        micro1 = micro1[:room]
      T.extend(micro1)

  # -------------------------
  # Micro-phase 2: alternate windows with pin-shrinking (guarded)
  # -------------------------
  if enable_alt_microphase and len(T) < CAP - 8:
    glo2 = min(l for l, r in T)
    ghi2 = max(r for l, r in T)
    G2 = max(1, ghi2 - glo2)

    seed_sz2 = max(8, min(24, len(T) // 400))
    U2 = _thin_seed(T, seed_sz2)
    if U2:
      ulo2 = min(l for l, r in U2)
      windows2 = [(0.06, 0.14), (0.28, 0.38), (0.54, 0.64), (0.78, 0.88)]
      blocks2 = []
      for (fa, fb) in windows2:
        win_lo2 = glo2 + int(round(fa * G2))
        base2 = win_lo2 - ulo2
        # Shrink to short "pins" to avoid omega spikes
        # Construct pins around midpoints with tiny span
        eps = max(1, G2 // 600)
        pins = []
        for idx, (l, r) in enumerate(U2):
          mid = int((l + r) // 2)
          L = mid + base2 + (idx % 3)  # stagger
          R = L + eps
          pins.append((L, R))
        blocks2.append(pins)

      # Distinct interleaving order
      micro2 = _interleave_blocks(blocks2, reverse=True, rotate=3)

      # Short pin connectors to tie windows
      eps2 = max(1, G2 // 1024)
      for pf in (0.10, 0.33, 0.66, 0.87):
        a = glo2 + int(round(pf * G2))
        micro2.append((a, a + eps2))

      room2 = CAP - len(T)
      if room2 > 0:
        if len(micro2) > room2:
          micro2 = micro2[:room2]
        T.extend(micro2)

  # -------------------------
  # Long-range connectors (cross-scale, conservative)
  # -------------------------
  if len(T) < CAP - 2:
    lo, hi, _ = _span(T)
    lr = _long_range_connectors(lo, hi)
    room3 = CAP - len(T)
    if room3 > 0 and lr:
      if len(lr) > room3:
        lr = lr[:room3]
      T.extend(lr)

  # Normalize to integer open intervals
  return _normalize(T)

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()