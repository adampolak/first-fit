# EVOLVE-BLOCK-START

def _span_delta(T):
  """Return (lo, hi, delta) for current interval set T."""
  lo = min(l for l, _ in T)
  hi = max(r for _, r in T)
  delta = hi - lo
  if delta <= 0:
    delta = 1
  return lo, hi, delta

def _build_blocks_sequential(T, starts, delta, lo, K=1):
  """
  Build four translated blocks sequentially (no interleaving).
  This ordering is designed to pressure FF without inflating omega.
  """
  S = []
  for s in starts:
    base = s * delta * K - lo
    S.extend((l + base, r + base) for (l, r) in T)
  return S

def _append_connectors(S, starts, delta, add_cross4=False):
  """
  Append the classic four connectors; optionally add a single long-range cross4
  on the final round to improve pressure while keeping total size under CAP.
  """
  s0, s1, s2, s3 = starts
  S.append(((s0 - 1) * delta, (s1 - 1) * delta))  # left cap
  S.append(((s2 + 2) * delta, (s3 + 2) * delta))  # right cap
  S.append(((s0 + 2) * delta, (s2 - 1) * delta))  # cross 1
  S.append(((s1 + 2) * delta, (s3 - 1) * delta))  # cross 2
  if add_cross4:
    S.append(((s0 + 4) * delta, (s3 + 4) * delta))  # long-range cross (final-round only)

def _append_micro_tail(T, max_extra=8):
  """
  Append a tiny tail of long caps near the end of T. This boosts FF late with minimal omega growth.
  Caps are inserted near the end (rather than simply appended) to intersect many active colors.
  """
  if not T or max_extra <= 0:
    return T
  lo, hi, delta = _span_delta(T)
  d2 = max(1, delta // 4)
  micro = [
      (lo + 1 * d2, lo + 5 * d2),
      (lo + 3 * d2, lo + 8 * d2),
      (hi - 6 * d2, hi - 2 * d2),
      (hi - 8 * d2, hi - 3 * d2),
  ]
  for i, interval in enumerate(micro[:max_extra]):
    pos = len(T) - (i * 2 + 1)
    if pos < 0:
      T.append(interval)
    else:
      T.insert(pos, interval)
  return T

def construct_intervals(seed_count=1, phase2_iters=2, enable_parabolic=False):
  """
  Deterministic, modular construction of intervals to maximize FirstFit colors
  per unit of omega, with a clean structure and tight CAP guards.
  - seed_count: initial seed scaffolding (kept small to control omega).
  - phase2_iters: number of micro delta2 rounds (0..2).
  - enable_parabolic: optional lightweight parabolic micro-phase for extra FF pressure.
  Returns: list of (l, r) intervals (integers), in FirstFit presentation order.
  """
  CAP = 9800

  # Stage 1: backbone KT spine
  spine_starts = (2, 6, 10, 14)
  template_bank = [
    (2, 6, 10, 14),  # classic KT
    (1, 5, 9, 13),   # left-shifted
    (3, 7, 11, 15),  # right-shifted
    (4, 8, 12, 16),  # stretched-right
  ]

  # Seed
  if seed_count <= 1:
    T = [(0, 1)]
  else:
    step = 3
    T = [(i * step, i * step + 1) for i in range(seed_count)]
  # Main spine: 6 rounds maximum while respecting CAP
  def apply_round(current_T, starts, interleave=False, add_cross4=False, rev_order=False):
    lo, hi, delta = _span_delta(current_T)
    blocks = []
    for s in starts:
      base = s * delta - lo
      blocks.append([(l + base, r + base) for (l, r) in current_T])

    S = []
    if interleave:
      maxlen = max((len(b) for b in blocks), default=0)
      order = list(range(len(blocks)))
      if rev_order:
        order.reverse()
      for i in range(maxlen):
        for idx in order:
          blk = blocks[idx]
          if i < len(blk):
            S.append(blk[i])
    else:
      if rev_order:
        blocks = list(reversed(blocks))
      for blk in blocks:
        S.extend(blk)

    _append_connectors(S, starts, delta, add_cross4=add_cross4)
    return S

  depth = 6
  for ridx in range(depth):
    starts = template_bank[ridx % len(template_bank)]
    interleave = (ridx % 2 == 0)  # even rounds interleave
    rev_order = False  # keep simple ordering pattern
    T = apply_round(T, starts, interleave=interleave, rev_order=rev_order)
    if len(T) >= CAP - 12:
      return T

  # Stage 2: micro delta2 rounds
  def thin_seed(current_T, max_seed):
    n = len(current_T)
    if n == 0 or max_seed <= 0:
      return []
    step = max(1, n // max_seed)
    return current_T[::step][:max_seed]

  def micro_round(current_T, round_id, budget):
    if budget <= 0 or not current_T:
      return []

    lo = min(l for l, _ in current_T)
    hi = max(r for _, r in current_T)
    G = max(1, hi - lo)
    delta2 = max(1, G // 2)

    per_block_target = max(8, min(64, budget // 12))
    U = thin_seed(current_T, per_block_target)
    if not U:
      return []

    ulo = min(l for l, _ in U)
    # Four windows (A realistic spread)
    window_fracs = [(0.12, 0.22), (0.35, 0.45), (0.58, 0.68), (0.80, 0.90)]
    blocks = []
    for (fa, fb) in window_fracs:
      win_lo = lo + int(round(fa * G))
      base = win_lo - ulo
      blocks.append([(l + base, r + base) for (l, r) in U])

    micro = []
    maxlen = max((len(b) for b in blocks), default=0)
    for i in range(maxlen):
      for blk in blocks:
        if i < len(blk):
          micro.append(blk[i])

    micro_connectors = [
      (lo + int(round(0.08 * G)), lo + int(round(0.30 * G))),
      (lo + int(round(0.60 * G)), lo + int(round(0.92 * G))),
      (lo + int(round(0.26 * G)), lo + int(round(0.56 * G))),
      (lo + int(round(0.44 * G)), lo + int(round(0.78 * G))),
    ]
    for a, b in micro_connectors:
      if b > a:
        micro.append((a, b))

    if len(micro) > budget:
      micro = micro[:budget]
    return micro

  remaining = CAP - len(T)
  steps = min(max(0, int(phase2_iters)), 2)
  for iter_id in range(steps):
    budget = CAP - len(T)
    if budget <= 0:
      break
    micro = micro_round(T, round_id=iter_id, budget=budget)
    if not micro:
      break
    avail = CAP - len(T)
    if len(micro) > avail:
      micro = micro[:avail]
    T.extend(micro)

  # Optional parabolic micro-phase
  if enable_parabolic:
    room = CAP - len(T)
    if room > 0:
      lo = min(l for l, _ in T)
      hi = max(r for _, r in T)
      delta = hi - lo if hi > lo else 1
      cap1 = (lo + max(1, int(round(0.15 * delta))), lo + max(1, int(round(0.60 * delta))))
      cap2 = (lo + max(1, int(round(0.25 * delta))), lo + max(1, int(round(0.80 * delta))))
      cap3 = (lo + max(1, int(round(0.55 * delta))), lo + max(1, int(round(0.95 * delta))))
      parabolic = [cap for cap in (cap1, cap2, cap3) if cap[1] > cap[0]]
      if parabolic:
        if len(parabolic) > room:
          parabolic = parabolic[:room]
        T.extend(parabolic)

  return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()