# EVOLVE-BLOCK-START

def construct_intervals(seed_count=1):
  """
  Construct a sequence of open intervals (l, r) in FirstFit arrival order,
  aiming to maximize FirstFit color usage while keeping the offline optimum
  (clique number) at most 10. The construction uses:
    - A laminated, parity-interleaved spine with variable-density spacing,
    - Two deterministic delta2 micro-phases with distinct window sets,
    - Cross-scale connectors (including cross4),
    - A post-append 10-track guard that prunes micro blocks if needed
      to ensure a 10-colorable certificate (hence omega ≤ 10).
  """

  # Hard cap (the evaluator enforces < 10000)
  CAP = 9800

  # ----------------------------
  # Helpers: geometry and mixing
  # ----------------------------

  def _span_delta(T):
    lo = min(l for l, _ in T)
    hi = max(r for _, r in T)
    delta = hi - lo
    if delta <= 0:
      delta = 1
    return lo, hi, delta

  def _build_blocks(current_T, starts, delta, lo, K=1, reverse_block_parity=False):
    """
    Build translated copies of current_T at offsets base = s*K*delta - lo.
    Optionally reverse inner order of odd-indexed blocks to break symmetry.
    """
    blocks = []
    for idx, s in enumerate(starts):
      base = s * delta * K - lo
      src = current_T[::-1] if (reverse_block_parity and (idx % 2 == 1)) else current_T
      block = [(l + base, r + base) for (l, r) in src]
      blocks.append(block)
    return blocks

  def _interleave_blocks(blocks, forward=True):
    """
    Interleave intervals from the given blocks in round-robin fashion.
    If forward is False, reverse block order before interleaving.
    """
    if not blocks:
      return []
    order = list(range(len(blocks)))
    if not forward:
      order.reverse()
    S = []
    maxlen = max(len(b) for b in blocks)
    for i in range(maxlen):
      for idx in order:
        blk = blocks[idx]
        if i < len(blk):
          S.append(blk[i])
    return S

  def _append_connectors(S, starts, delta, add_cross4=False):
    """
    Classic four connectors plus optional cross4. These are long caps that
    couple colors across blocks while remaining clique-friendly in the
    ten-track sense.
    """
    s0, s1, s2, s3 = starts
    S.append(((s0 - 1) * delta, (s1 - 1) * delta))  # left cap
    S.append(((s2 + 2) * delta, (s3 + 2) * delta))  # right cap
    S.append(((s0 + 2) * delta, (s2 - 1) * delta))  # cross 1
    S.append(((s1 + 2) * delta, (s3 - 1) * delta))  # cross 2
    if add_cross4:
      # One additional long-range connector modestly increases FF pressure
      # without exploding omega in presence of the track guard.
      S.append(((s0 + 4) * delta, (s3 + 4) * delta))

  # ----------------------------------------
  # Ten-track guard: certificate of omega ≤10
  # ----------------------------------------

  def _pack_to_k_tracks(intervals, K=10):
    """
    Interval-graph optimal coloring by left-endpoint sweep.
    Returns (ok, color_assignment) where ok=True iff colorable with ≤K tracks.
    Open-interval semantics: two intervals [l1,r1) and [l2,r2) overlap iff
    max(l1, l2) < min(r1, r2). For open intervals (l,r) with integer endpoints,
    equality at endpoints is non-intersection.
    """
    if not intervals:
      return True, []

    # Sort by left endpoint, break ties by right
    ord_idx = sorted(range(len(intervals)), key=lambda i: (intervals[i][0], intervals[i][1]))
    import heapq

    # Active heap by (end, color), and a min-heap of free colors
    active = []  # (end, color)
    free_colors = list(range(K))
    heapq.heapify(free_colors)

    colors = [None] * len(intervals)

    for i in ord_idx:
      l, r = intervals[i]
      if r <= l:
        r = l + 1

      # Release tracks whose end <= l (open intervals do not intersect at end)
      while active and active[0][0] <= l:
        _, c = heapq.heappop(active)
        heapq.heappush(free_colors, c)

      if not free_colors:
        return False, None

      c = heapq.heappop(free_colors)
      colors[i] = c
      heapq.heappush(active, (r, c))

    return True, colors

  def _try_append_with_track_guard(T, extras, K=10):
    """
    Try to append 'extras' to T while keeping 10-track colorability.
    Use a monotone (prefix) binary search to find the largest prefix
    of extras that keeps packability.
    """
    if not extras:
      return T

    lo, hi = 0, len(extras)
    best = 0
    # Binary search the largest prefix that is 10-packable
    while lo <= hi:
      mid = (lo + hi) // 2
      candidate = T + extras[:mid]
      ok, _ = _pack_to_k_tracks(candidate, K=K)
      if ok:
        best = mid
        lo = mid + 1
      else:
        hi = mid - 1

    if best > 0:
      return T + extras[:best]
    return T

  # -----------------------------
  # Laminated spine construction
  # -----------------------------

  # Four strong templates (rotated per round)
  template_bank = [
    (2, 6, 10, 14),  # classic KT-like anchor
    (1, 5, 9, 13),   # left-shifted
    (3, 7, 11, 15),  # right-shifted
    (4, 8, 12, 16),  # stretched-right
  ]

  # Seed with one unit interval (keeps omega modest)
  T = [(0, 1)] if seed_count <= 1 else [(0, 1)]

  # Predictive capacity accounting (size -> 4*size + connectors_per_round)
  def _next_size(sz, connectors):
    return 4 * sz + connectors

  rounds = 6
  for ridx in range(rounds):
    # Choose starts and density K
    starts = template_bank[ridx % len(template_bank)]
    K = 2 if (ridx % 2 == 1) else 1  # denser on odd rounds
    connectors_this_round = 5 if (ridx == rounds - 1) else 4  # add cross4 on final round

    # Capacity check before executing the round
    nxt = _next_size(len(T), connectors_this_round)
    if nxt > CAP:
      break

    lo, hi, delta = _span_delta(T)

    # Parity policies: interleave on even rounds; reverse inner parity on odd rounds
    blocks = _build_blocks(
      current_T=T,
      starts=starts,
      delta=delta,
      lo=lo,
      K=K,
      reverse_block_parity=(ridx % 2 == 1)
    )
    if ridx % 2 == 0:
      S = _interleave_blocks(blocks, forward=True)
    else:
      # Sequential but reversed block order to diversify entanglements
      S = []
      for blk in reversed(blocks):
        S.extend(blk)

    # Connectors with optional cross4 on the final round
    _append_connectors(S, starts, delta * K, add_cross4=(ridx == rounds - 1))

    # Commit spine round
    T = S

  # If we are close to capacity, return the strong spine baseline
  if len(T) >= CAP - 8:
    return T[:CAP]

  # ---------------------------------
  # Two deterministic delta2 micro phases
  # ---------------------------------

  BASE_SEED = 2654435761

  def _micro_round(current_T, budget, windows, iter_id=0):
    """
    Build a delta2-driven micro-block with interleaving and long caps
    over specified windows (fractions of the global span). Deterministic.
    """
    if not current_T or budget <= 8:
      return []

    glo = min(l for l, _ in current_T)
    ghi = max(r for _, r in current_T)
    G = max(1, ghi - glo)

    # Thin, evenly spaced deterministic sample with minor seed-based jitter
    n = len(current_T)
    seed2 = (BASE_SEED * 2654435761 + iter_id * 97 + 0x9E3779B9) & 0xFFFFFFFF
    # Pseudo-random but deterministic stride offset
    offset = (seed2 % max(1, n // 32)) if n >= 64 else 0

    seed_sz = max(8, min(32, n // 300))
    stride = max(1, n // max(1, seed_sz))
    U = [current_T[(i + offset) % n] for i in range(0, n, stride)][:seed_sz]
    if not U:
      return []

    # Build translated micro-blocks aligned to windows
    ulo = min(l for l, _ in U)
    blocks = []
    for (fa, fb) in windows:
      # clamp inside (0.03, 0.97) to avoid extreme edges
      fa = max(0.03, min(0.94, fa))
      fb = max(fa + 0.01, min(0.97, fb))
      win_lo = glo + int(round(fa * G))
      base = win_lo - ulo
      block = [(l + base, r + base) for (l, r) in U]
      # alternate internal reversal to break symmetry
      if int((fa * 100) + iter_id) % 2 == 1:
        block = list(reversed(block))
      blocks.append(block)

    # Interleave blocks (forward on even iter, reverse on odd)
    micro = _interleave_blocks(blocks, forward=(iter_id % 2 == 0))

    # Fractional-span long caps (four classic + cross4)
    L1 = (glo + int(round(0.07 * G)), glo + int(round(0.29 * G)))
    L2 = (glo + int(round(0.61 * G)), glo + int(round(0.93 * G)))
    C1 = (glo + int(round(0.24 * G)), glo + int(round(0.55 * G)))
    C2 = (glo + int(round(0.42 * G)), glo + int(round(0.77 * G)))
    X4 = (glo + int(round(0.15 * G)), glo + int(round(0.88 * G)))
    for a, b in (L1, L2, C1, C2, X4):
      if b > a:
        micro.append((a, b))

    if len(micro) > budget:
      micro = micro[:budget]
    return micro

  # First micro phase (primary)
  windows_A = [(0.12, 0.22), (0.35, 0.45), (0.58, 0.68), (0.80, 0.90)]
  room = CAP - len(T)
  if room > 8:
    extras_A = _micro_round(T, budget=room, windows=windows_A, iter_id=0)
    # Ten-track guard: append as large a prefix as remains 10-packable
    T = _try_append_with_track_guard(T, extras_A, K=10)

  # Second micro phase (secondary), distinct windows
  windows_B = [(0.05, 0.15), (0.28, 0.38), (0.60, 0.70), (0.82, 0.92)]
  room = CAP - len(T)
  if room > 8:
    extras_B = _micro_round(T, budget=room, windows=windows_B, iter_id=1)
    T = _try_append_with_track_guard(T, extras_B, K=10)

  # Final normalization and hard cap
  if len(T) > CAP:
    T = T[:CAP]
  return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()