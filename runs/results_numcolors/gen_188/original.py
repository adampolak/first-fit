# EVOLVE-BLOCK-START

def _span_delta(T):
  lo = min(l for l, _ in T)
  hi = max(r for _, r in T)
  delta = hi - lo
  if delta <= 0:
    delta = 1
  return lo, hi, delta

def _build_blocks_sequential(T, starts, delta, lo, K=1):
  """
  Build four translated blocks sequentially. No interleaving. This ordering
  is known to amplify FirstFit color pressure. K is a tunable spacing multiplier (default 1).
  """
  S = []
  for s in starts:
    base = s * delta * K - lo
    S.extend((l + base, r + base) for (l, r) in T)
  return S

def _append_connectors(S, starts, delta, add_cross4=False):
  """
  Append the classic four connectors. Optionally add a single long-range cross4
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
  Append a tiny tail of long caps near the end of T. This boosts FF late with minimal
  impact on omega. Caps are inserted near the end (rather than simply appended) to
  intersect many active colors while avoiding a single dense core.
  """
  if not T or max_extra <= 0:
    return T
  lo, hi, delta = _span_delta(T)
  d2 = max(1, delta // 4)
  # Symmetric, sparse long caps
  micro = [
      (lo + 1 * d2, lo + 5 * d2),
      (lo + 3 * d2, lo + 8 * d2),
      (hi - 6 * d2, hi - 2 * d2),
      (hi - 8 * d2, hi - 3 * d2),
  ]
  # Insert near the tail in a staggered fashion
  for i, interval in enumerate(micro[:max_extra]):
    pos = len(T) - (i * 2 + 1)
    if pos < 0:
      T.append(interval)
    else:
      T.insert(pos, interval)
  return T

def construct_intervals(seed_count=1, phase2_iters=2, enable_parabolic=False):
  """
  Construct a sequence of intervals of real line,
  in the order in which they are presented to FirstFit,
  so that it maximizes the number of colors used by FirstFit
  divided by the maximum number of intervals that cover a single point.

  This hybrid combines the six-round KT spine with delta2 micro-rounds
  and an optional parabolic micro-phase for extra FF pressure under control.
  """
  CAP = 9800

  # Stage 1: six KT rounds, deterministic spine with classic connectors
  spine_starts = (2, 6, 10, 14)
  T = [(0, 1)]

  def apply_round(current_T, starts, do_interleave=False):
    lo, hi, delta = _span_delta(current_T)
    blocks = []
    for s in starts:
      base = s * delta - lo
      blocks.append([(l + base, r + base) for (l, r) in current_T])

    S = []
    if do_interleave:
      maxlen = max(len(b) for b in blocks)
      for i in range(maxlen):
        for blk in blocks:
          if i < len(blk):
            S.append(blk[i])
    else:
      for blk in blocks:
        S.extend(blk)

    # connectors
    _append_connectors(S, starts, delta, add_cross4=False)
    return S

  for ridx in range(6):
    nxt_size = 4 * len(T) + 4
    if nxt_size > CAP:
      break
    # alternate interleaving pattern to push FF
    T = apply_round(T, spine_starts, do_interleave=(ridx % 2 == 0))

  # If we are near cap, return baseline
  if len(T) >= CAP - 16:
    return T

  # Stage 2: micro delta2 rounds with a thin seed
  def thin_seed(current_T, max_seed):
    n = len(current_T)
    if n == 0 or max_seed <= 0:
      return []
    step = max(1, n // max_seed)
    U = current_T[::step][:max_seed]
    return U

  def micro_round(current_T, round_id, budget):
    if budget <= 0 or not current_T:
      return []

    glo = min(l for l, r in current_T)
    ghi = max(r for l, r in current_T)
    G = max(1, ghi - glo)
    delta2 = max(1, G // 2)

    per_block_target = max(8, min(64, budget // 12))
    U = thin_seed(current_T, per_block_target)
    if not U:
      return []

    starts = spine_starts
    blocks = []
    ulo = min(l for l, r in U)
    for s in starts:
      base = s * delta2 - ulo
      block = [(l + base, r + base) for (l, r) in U]
      # tiny asymmetry
      if ((s // 2) % 2) == (round_id % 2):
        block = list(reversed(block))
      blocks.append(block)

    micro = []
    if round_id % 2 == 0:
      maxlen = max(len(b) for b in blocks)
      for i in range(maxlen):
        for blk in blocks:
          if i < len(blk):
            micro.append(blk[i])
    else:
      maxlen = max(len(b) for b in blocks)
      for i in range(maxlen):
        for blk in reversed(blocks):
          if i < len(blk):
            micro.append(blk[i])

    # connectors at delta2 scale
    s0, s1, s2, s3 = spine_starts
    connectors = [
      ((s0 - 1) * delta2 + glo, (s1 - 1) * delta2 + glo),
      ((s2 + 2) * delta2 + glo, (s3 + 2) * delta2 + glo),
      ((s0 + 2) * delta2 + glo, (s2 - 1) * delta2 + glo),
      ((s1 + 2) * delta2 + glo, (s3 - 1) * delta2 + glo),
      ((s0 + 3) * delta2 + glo, (s3 + 3) * delta2 + glo),
    ]
    micro.extend(connectors)

    cap1 = (glo + (delta2 // 4), glo + int(1.8 * delta2))
    cap2 = (glo + int(0.9 * delta2), glo + int(2.6 * delta2))
    mid = glo + G // 2
    cap3 = (mid - max(1, delta2 // 8), mid + max(1, delta2 // 8))
    for cap in (cap1, cap2, cap3):
      if cap[1] > cap[0]:
        micro.append(cap)

    if len(micro) > budget:
      micro = micro[:budget]
    return micro

  remaining = CAP - len(T)
  mr1 = micro_round(T, round_id=0, budget=max(0, remaining // 2))
  if mr1:
    room = CAP - len(T)
    if len(mr1) > room:
      mr1 = mr1[:room]
    T.extend(mr1)

  remaining = CAP - len(T)
  mr2 = micro_round(T, round_id=1, budget=max(0, remaining))
  if mr2:
    room = CAP - len(T)
    if len(mr2) > room:
      mr2 = mr2[:room]
    T.extend(mr2)

  # Optional parabolic micro-phase for extra diversity (lightweight)
  if enable_parabolic:
    room = CAP - len(T)
    if room > 0:
      lo = min(l for l, r in T)
      hi = max(r for l, r in T)
      delta = hi - lo if hi > lo else 1
      cap1 = (lo + max(1, int(round(0.15 * delta))), lo + max(1, int(round(0.60 * delta))))
      cap2 = (lo + max(1, int(round(0.25 * delta))), lo + max(1, int(round(0.80 * delta))))
      cap3 = (lo + max(1, int(round(0.55 * delta))), lo + max(1, int(round(0.95 * delta))))
      parabolic = [cap for cap in (cap1, cap2, cap3) if cap[1] > cap[0]]
      if parabolic:
        if len(parabolic) > room:
          parabolic = parabolic[:room]
        T.extend(parabolic)

  # Final return
  return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()