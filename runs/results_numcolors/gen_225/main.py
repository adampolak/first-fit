# EVOLVE-BLOCK-START

def construct_intervals(seed_count=1):
  """
  Construct a sequence of open intervals for FirstFit to force many colors
  while keeping the clique number <= 10. Returns a list of (l, r) integer pairs.

  Inputs:
    seed_count (int): retained for interface compatibility (unused beyond single-seed spine).
  Outputs:
    intervals: list[(int, int)]
  """

  # -----------------------------
  # Tunable parameters (safe defaults)
  # -----------------------------
  CAP = 9800               # hard capacity guard
  ROUNDS = 6               # KT-like main spine depth
  K_DENSITY = 2            # {1,2} densifies with thin "ghost" couplings when 2
  BASE_SEED = 137213       # deterministic seed for window jitter and parity
  CROSS4_ENABLED = True    # enable long-range cross connectors
  CROSS4_ROUND_PERIOD = 3  # per-round cross4 injection frequency
  MICRO_ITERS_A = 1        # first micro-phase iterations
  MICRO_ITERS_B = 1        # second micro-phase iterations

  # Window families for dual micro phases
  WINDOW_SET_A = [(0.12, 0.22), (0.35, 0.45), (0.58, 0.68), (0.80, 0.90)]
  WINDOW_SET_B = [(0.05, 0.15), (0.28, 0.38), (0.60, 0.70), (0.82, 0.92)]

  # -----------------------------
  # Small deterministic RNG (LCG)
  # -----------------------------
  def lcg_next(x):
    return (1103515245 * x + 12345) & 0x7fffffff

  def rnd01(x):
    x = lcg_next(x)
    return x, (x / 0x7fffffff)

  def jitter_pair(seed, a, b, scale=0.01):
    # Apply a tiny inward jitter so endpoints stay valid and inside span
    seed, ra = rnd01(seed)
    seed, rb = rnd01(seed)
    da = (ra - 0.5) * scale
    db = (rb - 0.5) * scale
    # Keep order a < b and remain within sensible interior
    aj = max(0.03, min(0.92, a + da))
    bj = max(0.08, min(0.97, b + db))
    if bj <= aj:
      bj = min(0.98, aj + 0.05)
    return seed, (aj, bj)

  # -----------------------------
  # Utility helpers
  # -----------------------------
  def _span(T):
    lo = min(l for l, r in T)
    hi = max(r for l, r in T)
    delta = hi - lo
    if delta <= 0:
      delta = 1
    return lo, hi, delta

  def _apply_connectors(S, starts, delta):
    s0, s1, s2, s3 = starts
    # Classic 4 connectors: left cap, right cap, two crosses
    S.append(((s0 - 1) * delta, (s1 - 1) * delta))
    S.append(((s2 + 2) * delta, (s3 + 2) * delta))
    S.append(((s0 + 2) * delta, (s2 - 1) * delta))
    S.append(((s1 + 2) * delta, (s3 - 1) * delta))

  def _apply_cross4(S, starts, delta):
    # Long-range connector tying distant blocks (~4 deltas apart)
    s0, s1, s2, s3 = starts
    S.append(((s0 + 4) * delta, (s3 + 4) * delta))

  def _thin_seed(current_T, max_seed):
    n = len(current_T)
    if n == 0 or max_seed <= 0:
      return []
    step = max(1, n // max_seed)
    return current_T[::step][:max_seed]

  # -----------------------------
  # Back bone templates
  # -----------------------------
  template_bank = [
    (2, 6, 10, 14),  # classic KT
    (1, 5, 9, 13),   # left-shifted
    (3, 7, 11, 15),  # right-shifted
    (4, 8, 12, 16),  # stretched-right
  ]

  # -----------------------------
  # Start with a single unit interval
  # -----------------------------
  T = [(0, 1)] if seed_count <= 1 else [(0, 1)]

  # -----------------------------
  # Round application with parity interleaving and K ghost density
  # -----------------------------
  def _apply_round(current_T, starts, round_id, k_density=1, cross4=False, parity_seed=0):
    lo, hi, delta = _span(current_T)

    # Build 4 translated full blocks
    blocks = []
    for s in starts:
      base = s * delta - lo
      block = [(l + base, r + base) for (l, r) in current_T]
      blocks.append(block)

    # Optional thin ghost couplings (built from a small seed to avoid heavy growth)
    ghost_blocks = []
    if k_density >= 2:
      seed_sz = max(8, min(48, len(current_T) // 240))
      U = _thin_seed(current_T, seed_sz)
      if U:
        ulo = min(l for l, r in U)
        # Place ghosts between consecutive main starts to couple colors at subscale
        for idx in range(len(starts) - 1):
          mid = (starts[idx] + starts[idx + 1]) // 2
          base = mid * delta - ulo
          gblock = [(l + base, r + base) for (l, r) in U]
          # Reverse every other ghost block to mix orderings
          if (idx + round_id + parity_seed) % 2 == 1:
            gblock = list(reversed(gblock))
          ghost_blocks.append(gblock)

    # Assemble with parity interleaving: even rounds forward, odd rounds reverse
    S = []
    all_blocks = list(blocks)  # keep main blocks leading
    # Insert ghosts interleaved across blocks for mixing, but in a bounded way
    interleave_sets = all_blocks + ghost_blocks
    forward = (round_id % 2 == 0)
    if not forward:
      interleave_sets = list(reversed(interleave_sets))
    maxlen = max(len(b) for b in interleave_sets)
    for i in range(maxlen):
      for blk in interleave_sets:
        if i < len(blk):
          S.append(blk[i])

    # Add connectors, classic + optional cross4
    _apply_connectors(S, starts, delta)
    if cross4:
      _apply_cross4(S, starts, delta)

    return S

  # -----------------------------
  # Stage 1: KT-like spine with deterministic parity and K=2 density
  # -----------------------------
  # Predictive size guard: sz -> 4*sz + 4 (+ ghosts, connectors)
  for ridx in range(ROUNDS):
    nxt = 4 * len(T) + 4
    if nxt > CAP:
      break
    starts = template_bank[ridx % len(template_bank)]
    # Cross4 enabled every CROSS4_ROUND_PERIOD-th round (well below CAP)
    add_cross4 = CROSS4_ENABLED and (ridx % CROSS4_ROUND_PERIOD == 1)
    T = _apply_round(
      T,
      starts,
      round_id=ridx,
      k_density=K_DENSITY,
      cross4=add_cross4,
      parity_seed=(BASE_SEED ^ (ridx * 1315423911)) & 0x7fffffff
    )
    if len(T) >= CAP:
      T = T[:CAP]
      return T

  # Early exit if near capacity
  if len(T) >= CAP - 16:
    return T

  # -----------------------------
  # Micro-phase builder (parametric)
  # -----------------------------
  def _build_micro(current_T, budget, window_set, iter_id, base_seed, add_long_cross=False):
    if not current_T or budget <= 8:
      return []

    glo = min(l for l, r in current_T)
    ghi = max(r for l, r in current_T)
    G = max(1, ghi - glo)

    # Seed size is thin to avoid raising omega; slight growth with iteration id
    seed_goal = max(10, min(56, len(current_T) // (260 - 10 * (iter_id % 3))))
    U = _thin_seed(current_T, seed_goal)
    if not U:
      return []

    ulo = min(l for l, r in U)

    # Apply tiny deterministic jitter to each window
    seed = (base_seed ^ (iter_id * 2654435761)) & 0x7fffffff
    windows = []
    for (a, b) in window_set:
      seed, (aj, bj) = jitter_pair(seed, a, b, scale=0.012)
      windows.append((aj, bj))

    # Build translated micro-blocks aligned to these windows
    blocks = []
    for idx, (fa, fb) in enumerate(windows):
      win_lo = glo + int(round(fa * G))
      base = win_lo - ulo
      block = [(l + base, r + base) for (l, r) in U]
      # Internal reversal pattern to diversify orderings deterministically
      if ((idx + iter_id) ^ (base_seed & 3)) % 2 == 1:
        block = list(reversed(block))
      blocks.append(block)

    # Interleave with parity dependent on iter_id and seed
    order = list(range(len(blocks)))
    if ((iter_id + (base_seed & 7)) % 2) == 1:
      order.reverse()
    micro = []
    maxlen = max(len(blocks[i]) for i in order)
    for i in range(maxlen):
      for j in order:
        blk = blocks[j]
        if i < len(blk):
          micro.append(blk[i])

    # Fractional-span connectors (scaled to G)
    micro_connectors = [
      (glo + int(round(0.08 * G)), glo + int(round(0.30 * G))),
      (glo + int(round(0.60 * G)), glo + int(round(0.92 * G))),
      (glo + int(round(0.26 * G)), glo + int(round(0.56 * G))),
      (glo + int(round(0.44 * G)), glo + int(round(0.78 * G))),
    ]
    if add_long_cross:
      micro_connectors.append((glo + int(round(0.18 * G)), glo + int(round(0.84 * G))))

    for a, b in micro_connectors:
      if b > a:
        micro.append((a, b))

    if len(micro) > budget:
      micro = micro[:budget]
    return micro

  # -----------------------------
  # Stage 2: Dual micro phases (A then B), CAP-guarded
  # -----------------------------
  # Micro A
  steps = max(0, int(MICRO_ITERS_A))
  for it in range(steps):
    room = CAP - len(T)
    if room <= 8:
      break
    microA = _build_micro(
      T, room, WINDOW_SET_A, iter_id=it,
      base_seed=(BASE_SEED ^ 0x9e3779b1),
      add_long_cross=False
    )
    if not microA:
      break
    if len(microA) > room:
      microA = microA[:room]
    T.extend(microA)

  if len(T) >= CAP - 16:
    return T

  # Micro B (alternate windows, reverse parity via seed, include long cross)
  stepsB = max(0, int(MICRO_ITERS_B))
  for it in range(stepsB):
    room = CAP - len(T)
    if room <= 8:
      break
    microB = _build_micro(
      T, room, WINDOW_SET_B, iter_id=it,
      base_seed=(BASE_SEED ^ 0x7f4a7c15),
      add_long_cross=True
    )
    if not microB:
      break
    if len(microB) > room:
      microB = microB[:room]
    T.extend(microB)

  # -----------------------------
  # Stage 3: Global cross4 connector layer (very sparse, CAP-gated)
  # -----------------------------
  if CROSS4_ENABLED and len(T) < CAP - 4:
    lo, hi, _ = _span(T)
    G = max(1, hi - lo)
    # Three sparse long connectors across distant span portions
    cross4_layer = [
      (lo + int(0.03 * G), lo + int(0.37 * G)),
      (lo + int(0.21 * G), lo + int(0.66 * G)),
      (lo + int(0.52 * G), lo + int(0.94 * G)),
    ]
    for a, b in cross4_layer:
      if len(T) >= CAP:
        break
      if b > a:
        T.append((a, b))

  # Final CAP trim
  if len(T) > CAP:
    T = T[:CAP]

  return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()