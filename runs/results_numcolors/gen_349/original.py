# EVOLVE-BLOCK-START

def construct_intervals(seed_count=1):
  """
  Deterministic fractal KT-style spine with a six-round backbone and two-stage micro-phases.
  Returns a list of open intervals (l, r) in FirstFit presentation order.
  """

  # Cap to keep total intervals under ~10000
  CAP = 9800

  # Backbone: four strong KT templates rotated across rounds
  template_bank = [
    (2, 6, 10, 14),  # classic KT
    (1, 5, 9, 13),   # left-shifted
    (3, 7, 11, 15),  # right-shifted
    (4, 8, 12, 16),  # stretched-right
  ]

  # Seed with a single unit interval to avoid early omega inflation
  T = [(0, 1)]

  # Helpers
  def _span_delta(current_T):
    lo = min(l for l, r in current_T)
    hi = max(r for l, r in current_T)
    delta = hi - lo
    if delta <= 0:
      delta = 1
    return lo, hi, delta

  def _append_connectors(S, starts, delta):
    # Classic four connectors that couple blocks without exploding omega
    s0, s1, s2, s3 = starts
    S.append(((s0 - 1) * delta, (s1 - 1) * delta))  # left cap
    S.append(((s2 + 2) * delta, (s3 + 2) * delta))  # right cap
    S.append(((s0 + 2) * delta, (s2 - 1) * delta))  # cross 1
    S.append(((s1 + 2) * delta, (s3 - 1) * delta))  # cross 2

  def _apply_round(current_T, starts, interleave=False, reverse_order=False, extra_connect=False):
    lo, hi, delta = _span_delta(current_T)

    # Build four translated blocks
    blocks = []
    for s in starts:
      base = s * delta - lo
      block = [(l + base, r + base) for (l, r) in current_T]
      blocks.append(block)

    # Assemble S with optional interleaving and reverse order
    S = []
    if interleave:
      maxlen = max(len(b) for b in blocks)
      order = list(range(4))
      if reverse_order:
        order.reverse()
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

    # Append connectors
    _append_connectors(S, starts, delta)

    # Optional extra cross (only when explicitly requested)
    if extra_connect:
      s0, s1, s2, s3 = starts
      S.append(((s0 + 4) * delta, (s3 + 4) * delta))

    return S

  # Stage 1: six KT rounds with deterministic rotation; interleaving boosts FF mixing
  depth = 6
  for ridx in range(depth):
    starts = template_bank[ridx % len(template_bank)]
    interleave = True
    rev = (ridx % 2 == 1)
    T = _apply_round(T, starts, interleave=interleave, reverse_order=rev)
    if len(T) >= CAP:
      T = T[:CAP]
      return T

  # Early exit guard
  if len(T) >= CAP - 8:
    return T

  # Stage 1b: tail-tail caps to tighten and broaden color coverage
  lo, hi, delta = _span_delta(T)
  span = max(1, hi - lo)
  def cap(a_frac, b_frac):
    L = lo + max(1, int(round(a_frac * span)))
    R = lo + max(1, int(round(b_frac * span)))
    if R <= L:
      R = L + 1
    return (L, R)
  tail_caps = [cap(0.08, 0.60), cap(0.25, 0.75), cap(0.75, 0.92)]
  room = CAP - len(T)
  if room > 0:
    for i, iv in enumerate(tail_caps[:room]):
      pos = len(T) - (i * 2 + 1)
      if pos < 0:
        T.append(iv)
      else:
        T.insert(pos, iv)

  if len(T) >= CAP - 16:
    return T

  # Stage 2: micro delta rounds (two families for robust mixing)
  def micro_round(current_T, budget, iter_id=0, alt=False):
    if not current_T or budget <= 8:
      return []

    lo = min(l for l, r in current_T)
    hi = max(r for l, r in current_T)
    delta = max(1, hi - lo)

    n = len(current_T)
    seed_sz = max(6, min(40, n // 200))
    stride = max(1, n // seed_sz)
    U = [current_T[i] for i in range(0, n, stride)][:seed_sz]
    if not U:
      return []

    ulo = min(l for l, r in U)

    # Two families of windows: primary (A) and alternate (B)
    if not alt:
      window_fracs = [(0.10, 0.20), (0.32, 0.42), (0.58, 0.68), (0.82, 0.92)]
    else:
      window_fracs = [(0.05, 0.15), (0.25, 0.35), (0.50, 0.60), (0.75, 0.85)]

    blocks = []
    for (fa, fb) in window_fracs:
      win_lo = lo + int(round(fa * delta))
      base = win_lo - ulo
      block = [(l + base, r + base) for (l, r) in U]
      blocks.append(block)

    micro = []
    maxlen = max(len(b) for b in blocks)
    order = list(range(len(blocks)))
    if iter_id % 2 == 1:
      order.reverse()
    for i in range(maxlen):
      for idx in order:
        blk = blocks[idx]
        if i < len(blk):
          micro.append(blk[i])

    # Cross-scale connectors (deterministic)
    connectors = [
      (lo + int(round(0.08 * delta)), lo + int(round(0.30 * delta))),
      (lo + int(round(0.60 * delta)), lo + int(round(0.92 * delta))),
      (lo + int(round(0.26 * delta)), lo + int(round(0.56 * delta))),
      (lo + int(round(0.44 * delta)), lo + int(round(0.78 * delta))),
    ]
    for a, b in connectors:
      if b > a:
        micro.append((a, b))

    if len(micro) > budget:
      micro = micro[:budget]
    return micro

  # Stage 2 micro-rounds (two rounds, primary then alternate)
  room = CAP - len(T)
  if room > 8:
    micro1 = micro_round(T, room, iter_id=0, alt=False)
    if micro1:
      if len(micro1) > room:
        micro1 = micro1[:room]
      T.extend(micro1)

  room = CAP - len(T)
  if room > 8:
    micro2 = micro_round(T, room, iter_id=1, alt=True)
    if micro2:
      if len(micro2) > room:
        micro2 = micro2[:room]
      T.extend(micro2)

  # Final cross-scale connectors anchored to current span
  room = CAP - len(T)
  if room > 0:
    lo2, hi2, delta2 = _span_delta(T)
    def frac_iv(a_frac, b_frac):
      L = lo2 + max(1, int(round(a_frac * delta2)))
      R = lo2 + max(1, int(round(b_frac * delta2)))
      if R <= L:
        R = L + 1
      return (L, R)
    connectors2 = [
      frac_iv(0.14, 0.52),
      frac_iv(0.32, 0.88),
      frac_iv(0.05, 0.47),
      frac_iv(0.62, 0.95),
      frac_iv(0.20, 0.41),
      frac_iv(0.71, 0.89),
    ]
    if room < len(connectors2):
      connectors2 = connectors2[:room]
    T.extend(connectors2)

  # Final cap trim
  if len(T) > CAP:
    T = T[:CAP]

  return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()