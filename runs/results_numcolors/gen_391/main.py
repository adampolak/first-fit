# EVOLVE-BLOCK-START

def construct_intervals(enable_alt_microphase=True):
  """
  Deterministic spine-wave with an expanded backbone (six templates),
  two micro-phases (A and optional B), and deterministic post-phase connectors.
  Keeps omega under a target bound while maximizing FirstFit color usage.
  """
  CAP = 9800

  # Expanded backbone bank of six templates
  TEMPLATE_BANK = [
    (2, 6, 10, 14),  # classic KT
    (1, 5, 9, 13),   # left-shifted
    (3, 7, 11, 15),  # right-shifted
    (4, 8, 12, 16),  # stretched-right
    (5, 9, 13, 17),  # extra-right
    (6, 10, 14, 18), # extra-classic
  ]

  # Seed spine with a single unit interval
  T = [(0, 1)]

  def span_params(seq):
    lo = min(l for l, r in seq)
    hi = max(r for l, r in seq)
    d = hi - lo
    if d <= 0:
      d = 1
    return lo, hi, d

  def apply_round(curr, starts, interleave=False, reverse_order=False):
    lo, hi, delta = span_params(curr)
    blocks = []
    for s in starts:
      base = s * delta - lo
      blocks.append([(l + base, r + base) for (l, r) in curr])

    S = []
    if interleave:
      order = list(range(4))
      if reverse_order:
        order = list(reversed(order))
      maxlen = max(len(b) for b in blocks)
      for i in range(maxlen):
        for idx in order:
          blk = blocks[idx]
          if i < len(blk):
            S.append(blk[i])
    else:
      if reverse_order:
        blocks = list(reversed(blocks))
      for blk in blocks:
        S.extend(blk)

    s0, s1, s2, s3 = starts
    S.append(((s0 - 1) * delta, (s1 - 1) * delta))  # left cap
    S.append(((s2 + 2) * delta, (s3 + 2) * delta))  # right cap
    S.append(((s0 + 2) * delta, (s2 - 1) * delta))  # cross 1
    S.append(((s1 + 2) * delta, (s3 - 1) * delta))  # cross 2
    return S

  # Stage 1: Backbone rounds
  for ridx in range(9):  # up to 9 rounds
    nxt = 4 * len(T) + 4
    if nxt > CAP:
      break
    starts = TEMPLATE_BANK[ridx % len(TEMPLATE_BANK)]
    interleave = (ridx % 2 == 0)
    reverse_order = (ridx % 2 == 1)
    T = apply_round(T, starts, interleave=interleave, reverse_order=reverse_order)
    if len(T) >= CAP - 16:
      break

  # Early exit if near CAP
  if len(T) >= CAP - 8:
    return normalize_and_cap(T)

  # Densification: safe splitting of long intervals to create more colors without increasing omega
  if len(T) < CAP - 40:
    candidates = [i for i, (l, r) in enumerate(T) if (r - l) > 6][:8]
    for idx in sorted(candidates, reverse=True):
      l, r = T[idx]
      mid = (l + r) // 2
      if mid > l:
        T[idx] = (l, mid)
        T.insert(idx + 1, (mid, r))
      if len(T) >= CAP:
        break

  lo, hi, span_ = span_params(T)

  def cap_at(a_frac, b_frac):
    L = lo + max(1, int(round(a_frac * span_)))
    R = lo + max(1, int(round(b_frac * span_)))
    if R <= L:
      R = L + 1
    return (L, R)

  # Micro-phase A: three long-range caps
  caps = [cap_at(0.08, 0.60), cap_at(0.25, 0.75), cap_at(0.75, 0.92)]
  for c in caps:
    if len(T) >= CAP:
      break
    if c[1] > c[0]:
      T.append(c)

  # Micro-phase builder for windowed micro-blocks
  def build_micro(curr, budget, windows, alt=False, long_cross=False, seed_boost=False):
    if not curr or budget <= 0:
      return []
    glo = min(l for l, r in curr)
    ghi = max(r for l, r in curr)
    G = max(1, ghi - glo)

    base_seed = max(6, min(40, len(curr) // 120))
    seed_sz = base_seed + (4 if seed_boost else 0)
    seed_sz = min(seed_sz, 60)

    stride = max(1, len(curr) // seed_sz)
    U = [curr[i] for i in range(0, len(curr), stride)][:seed_sz]
    if not U:
      return []

    ulo = min(l for l, r in U)

    blocks = []
    for (fa, fb) in windows:
      fa = max(0.05, min(0.95, fa))
      fb = max(0.10, min(0.95, fb))
      win_lo = glo + int(round(fa * G))
      base = win_lo - ulo
      block = [(l + base, r + base) for (l, r) in U]
      if alt:
        tag = (int(round(fa * 100)) // 5) % 2
        if tag == 1:
          block = list(reversed(block))
      blocks.append(block)

    micro = []
    maxlen = max((len(b) for b in blocks), default=0)
    if not blocks:
      return []
    if not alt:
      for i in range(maxlen):
        for blk in blocks:
          if i < len(blk):
            micro.append(blk[i])
    else:
      for i in range(maxlen):
        for blk in reversed(blocks):
          if i < len(blk):
            micro.append(blk[i])

    connectors = [
      (glo + int(round(0.08 * G)), glo + int(round(0.30 * G))),
      (glo + int(round(0.60 * G)), glo + int(round(0.92 * G))),
      (glo + int(round(0.26 * G)), glo + int(round(0.56 * G))),
      (glo + int(round(0.44 * G)), glo + int(round(0.78 * G))),
    ]
    for a, b in connectors:
      if b > a:
        micro.append((a, b))

    if long_cross:
      a = glo + int(round(0.18 * G))
      b = glo + int(round(0.84 * G))
      if b > a:
        micro.append((a, b))

    if len(micro) > budget:
      micro = micro[:budget]
    return micro

  room = CAP - len(T)
  if room > 8:
    budgetA = max(0, int(room * 0.60))
    microA = build_micro(T, budgetA, windows=[(0.12,0.22),(0.35,0.45),(0.58,0.68),(0.80,0.90)], alt=False, long_cross=False, seed_boost=False)
    if microA:
      avail = CAP - len(T)
      if len(microA) > avail: microA = microA[:avail]
      T.extend(microA)

  room = CAP - len(T)
  if enable_alt_microphase and room > 8:
    budgetB = room
    microB = build_micro(T, budgetB, windows=[(0.06,0.14),(0.28,0.38),(0.54,0.64),(0.74,0.84)], alt=True, long_cross=True, seed_boost=True)
    if microB:
      avail = CAP - len(T)
      if len(microB) > avail: microB = microB[:avail]
      T.extend(microB)

  # Post-micro deterministic connectors
  room = CAP - len(T)
  if room > 0:
    lr = []
    lo2, hi2, d2 = span_params(T)
    anchors = [
      (0.12, 0.62),
      (0.22, 0.78),
      (0.38, 0.86),
      (0.48, 0.90),
      (0.30, 0.70),
      (0.15, 0.85),
    ]
    rot = 0
    for i in range(min(room, len(anchors))):
      a, b = anchors[(rot + i) % len(anchors)]
      L = lo2 + max(1, int(round(a * d2)))
      R = lo2 + max(1, int(round(b * d2)))
      if R > L:
        lr.append((L, R))
    if lr:
      if len(lr) > room:
        lr = lr[:room]
      T.extend(lr)

  return normalize_and_cap(T)

def run_experiment(**kwargs):
  return construct_intervals(enable_alt_microphase=kwargs.get('enable_alt_microphase', True))

def normalize_and_cap(seq):
  if not seq:
    return []
  lo = min(l for l, r in seq)
  if lo < 0:
    seq = [(l - lo, r - lo) for (l, r) in seq]
  out = []
  for l, r in seq:
    li = int(round(l))
    ri = int(round(r))
    if ri <= li:
      ri = li + 1
    out.append((li, ri))
  if len(out) > 9800:
    out = out[:9800]
  return out

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()