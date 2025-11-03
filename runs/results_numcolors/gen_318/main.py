# EVOLVE-BLOCK-START

def construct_intervals(enable_alt_microphase=True):
  """
  Construct a sequence of open intervals for FirstFit to maximize FF/omega.
  Deterministic six-round spine with a seed-driven interleaving policy,
  expanded template bank, disciplined connectors, and two micro-phase passes.

  Args:
    enable_alt_microphase (bool): enable a second, alternate micro-phase pass.

  Returns:
    list[tuple[int, int]]: intervals (l, r) as open intervals, in arrival order.
  """

  # ------------------------------
  # Global configuration and guards
  # ------------------------------
  CAP = 9800               # hard capacity (< 10000)
  ROUNDS = 6               # spine depth
  BASE_SEED = 0x9E3779B1   # deterministic base seed (golden ratio-based)

  # Six four-start templates (expanded bank).
  TEMPLATE_BANK = [
    (2, 6, 10, 14),  # classic KT
    (1, 5, 9, 13),   # left-shifted
    (3, 7, 11, 15),  # right-shifted
    (4, 8, 12, 16),  # stretched-right
    (2, 5, 11, 14),  # stagger variant A
    (3, 6, 12, 15),  # stagger variant B
  ]

  # Seed: single unit interval
  T = [(0, 1)]

  # ------------------------------
  # Small utilities
  # ------------------------------
  def span(current_T):
    lo = min(l for l, r in current_T)
    hi = max(r for l, r in current_T)
    delta = hi - lo
    if delta <= 0:
      delta = 1
    return lo, hi, delta

  def lcg_step(s, inc):
    # Simple LCG mix; keep deterministic and platform-independent
    s = (s * 1664525 + 1013904223 + (inc & 0xFFFFFFFF)) & 0xFFFFFFFF
    return s

  def permute_order(seed, n=4):
    # Produce a permutation of range(n) from the seed with two swaps
    order = list(range(n))
    s = seed
    for k in range(2):
      s = lcg_step(s, 0xA5A5A5A5 ^ (k * 0x1F))
      i = (s >> 16) % n
      s = lcg_step(s, 0xC3EF1A55 ^ (k * 0x2B))
      j = (s >> 8) % n
      if i != j:
        order[i], order[j] = order[j], order[i]
    return order

  def interleave_blocks(blocks, order):
    # Round-robin interleaving in given block order
    out = []
    maxlen = max(len(b) for b in blocks)
    for i in range(maxlen):
      for idx in order:
        blk = blocks[idx]
        if i < len(blk):
          out.append(blk[i])
    return out

  def connectors_for_round(starts, delta, variant=0):
    # Several deterministic connector variants; all conservative for omega
    s0, s1, s2, s3 = starts
    conns = []
    if variant == 0:
      # Classic KT connectors (Figure 4 style)
      conns.append(((s0 - 1) * delta, (s1 - 1) * delta))
      conns.append(((s2 + 2) * delta, (s3 + 2) * delta))
      conns.append(((s0 + 2) * delta, (s2 - 1) * delta))
      conns.append(((s1 + 2) * delta, (s3 - 1) * delta))
    elif variant == 1:
      # Slightly widened caps; cross ties unchanged
      conns.append(((s0 - 2) * delta, (s1) * delta))
      conns.append(((s2 + 1) * delta, (s3 + 3) * delta))
      conns.append(((s0 + 2) * delta, (s2 - 1) * delta))
      conns.append(((s1 + 2) * delta, (s3 - 1) * delta))
    else:
      # Balanced, slightly shifted caps
      conns.append(((s0 - 1) * delta, (s1) * delta))
      conns.append(((s2 + 1) * delta, (s3 + 2) * delta))
      conns.append(((s0 + 2) * delta, (s2 - 1) * delta))
      conns.append(((s1 + 2) * delta, (s3 - 1) * delta))
    # Ensure endpoints monotone
    out = []
    for a, b in conns:
      if b <= a:
        b = a + 1
      out.append((a, b))
    return out

  # ------------------------------
  # Spine construction (six deterministic rounds)
  # ------------------------------
  seed = BASE_SEED
  for ridx in range(ROUNDS):
    # Predict capacity usage: size -> 4*size + 4 (connectors included later)
    if 4 * len(T) + 4 > CAP:
      break

    starts = TEMPLATE_BANK[ridx % len(TEMPLATE_BANK)]
    lo, hi, delta = span(T)

    # Build translated blocks for the four starts
    raw_blocks = []
    for s in starts:
      base = s * delta - lo
      block = [(l + base, r + base) for (l, r) in T]
      raw_blocks.append(block)

    # Deterministic block permutation and mandatory interleaving
    seed = lcg_step(seed, 0xDEADBEEF ^ ridx)
    order = permute_order(seed % (1 << 32), n=4)
    S = interleave_blocks(raw_blocks, order)

    # Deterministic connector variant per round
    conn_variant = (seed >> 29) % 3  # 0,1,2
    S.extend(connectors_for_round(starts, delta, variant=conn_variant))

    T = S

  # ------------------------------
  # Long-range connector discipline (post-spine)
  # ------------------------------
  if len(T) < CAP - 4:
    lo, hi, delta = span(T)
    G = max(1, hi - lo)
    # Deterministic fractional connectors across the global span
    lr_fracs = [
      (0.06, 0.94),
      (0.14, 0.86),
      (0.22, 0.78),
      (0.34, 0.66),
      (0.18, 0.84),
      (0.28, 0.72),
    ]
    LR = []
    for (fa, fb) in lr_fracs:
      a = lo + int(round(fa * G))
      b = lo + int(round(fb * G))
      if b <= a:
        b = a + 1
      LR.append((a, b))
    # Only append up to available capacity slack
    room = CAP - len(T)
    if room > 0:
      T.extend(LR[:room])

  # Early exit if near capacity
  if len(T) >= CAP - 8:
    return T if len(T) <= CAP else T[:CAP]

  # ------------------------------
  # Micro-phase builder (generic)
  # ------------------------------
  def build_micro(current_T, budget, windows, seed_tag):
    if not current_T or budget <= 8:
      return []

    glo = min(l for l, r in current_T)
    ghi = max(r for l, r in current_T)
    G = max(1, ghi - glo)

    # Thin, evenly spaced seed; deterministic size choice
    base_seed = max(8, min(32, len(current_T) // 300))
    seed_sz = base_seed
    stride = max(1, len(current_T) // max(1, seed_sz))
    U = [current_T[i] for i in range(0, len(current_T), stride)][:seed_sz]
    if not U:
      return []

    ulo = min(l for l, r in U)

    # Build translated micro-blocks aligned to windows
    blocks = []
    for (fa, fb) in windows:
      win_lo = glo + int(round(fa * G))
      base = win_lo - ulo
      block = [(l + base, r + base) for (l, r) in U]
      # Alternate internal reversal based on seed_tag and window index
      if int(round((fa + fb) * 100)) % 2 == (seed_tag % 2):
        block = list(reversed(block))
      blocks.append(block)

    # Deterministic interleaving order across blocks
    order = permute_order(lcg_step(BASE_SEED, seed_tag) % (1 << 32), n=len(blocks))
    # permute_order for n != 4: adapt by truncating/expanding
    def permute_custom(order, n):
      # order is over [0..3]; derive a permutation over [0..n-1]
      base = list(range(n))
      s = lcg_step(BASE_SEED ^ (seed_tag * 0x45D9F3B), n)
      for k in range(n + 1):
        i = (s >> (k % 17)) % n
        j = (s >> ((k + 7) % 19)) % n
        if i != j:
          base[i], base[j] = base[j], base[i]
        s = lcg_step(s, k * 0x9E37)
      return base

    block_order = permute_custom(order, len(blocks))
    micro = []
    maxlen = max(len(b) for b in blocks)
    for i in range(maxlen):
      for idx in block_order:
        blk = blocks[idx]
        if i < len(blk):
          micro.append(blk[i])

    # Fractional-span connectors within micro-phase
    micro_connectors = [
      (glo + int(round(0.08 * G)), glo + int(round(0.30 * G))),
      (glo + int(round(0.60 * G)), glo + int(round(0.92 * G))),
      (glo + int(round(0.26 * G)), glo + int(round(0.56 * G))),
      (glo + int(round(0.44 * G)), glo + int(round(0.78 * G))),
    ]
    for a, b in micro_connectors:
      if b > a:
        micro.append((a, b))

    # Trim to budget
    if len(micro) > budget:
      micro = micro[:budget]
    return micro

  # ------------------------------
  # Micro-phase Pass A (primary windows)
  # ------------------------------
  room = CAP - len(T)
  if room > 8:
    windows_A = [
      (0.12, 0.22),
      (0.35, 0.45),
      (0.58, 0.68),
      (0.80, 0.90),
    ]
    microA = build_micro(T, room, windows_A, seed_tag=1)
    if microA:
      T.extend(microA[: (CAP - len(T))])

  # ------------------------------
  # Micro-phase Pass B (independent windows)
  # ------------------------------
  room = CAP - len(T)
  if enable_alt_microphase and room > 8:
    windows_B = [
      (0.16, 0.26),
      (0.40, 0.50),
      (0.62, 0.72),
      (0.84, 0.94),
    ]
    microB = build_micro(T, room, windows_B, seed_tag=2)
    if microB:
      T.extend(microB[: (CAP - len(T))])

  # Final capacity trim
  if len(T) > CAP:
    T = T[:CAP]
  return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()