# EVOLVE-BLOCK-START

def construct_intervals():
  """
  Evolved KT-style spine construction with guarded micro-phases.

  Returns:
    intervals: list of (l, r) integer tuples (open intervals), in FF arrival order.
  """
  CAP = 9800  # hard limit for interval count

  # --- Configuration knobs (tunable) ---
  base_seed = 0xC0FFEE  # deterministic global seed
  use_density_K = 2     # try denser packing (1 or 2). K=2 often increases FF pressure.
  max_spine_rounds = 6  # target spine rounds (stopped earlier if capacity risk)
  # Four-start templates (rotated per round for diversity)
  template_bank = [
    (2, 6, 10, 14),
    (1, 5, 9, 13),
    (3, 7, 11, 15),
    (2, 4, 8, 12),
  ]

  # --- Utility helpers ---
  def span_delta(T):
    lo = min(l for l, _ in T)
    hi = max(r for _, r in T)
    d = hi - lo
    if d <= 0:
      d = 1
    return lo, hi, d

  def predict_next_size(sz):
    # KT-style recurrence: next = 4*sz + 4 connectors
    return 4 * sz + 4

  def max_rounds_within_cap(initial_size, max_rounds):
    sz = initial_size
    done = 0
    for _ in range(int(max_rounds)):
      nxt = predict_next_size(sz)
      if nxt > CAP:
        break
      sz = nxt
      done += 1
    return done, sz

  def lcg_next(x):
    # simple deterministic linear congruential step for seed derivation
    return (x * 2654435761 + 979) & 0x7FFFFFFF

  # --- Connectors utilities ---
  def append_classic_connectors(S, starts, delta, lo=0, add_cross4=False, scale=1):
    """Append classic 4 connectors (and optional cross4). scale multiplies the spacing."""
    s0, s1, s2, s3 = starts
    # localize to lo if requested: many micro connectors add glo offset already.
    S.append(((s0 - 1) * delta * scale + lo, (s1 - 1) * delta * scale + lo))
    S.append(((s2 + 2) * delta * scale + lo, (s3 + 2) * delta * scale + lo))
    S.append(((s0 + 2) * delta * scale + lo, (s2 - 1) * delta * scale + lo))
    S.append(((s1 + 2) * delta * scale + lo, (s3 - 1) * delta * scale + lo))
    if add_cross4:
      # long-range cross connector (keeps omega low if used sparingly)
      S.append(((s0 + 4) * delta * scale + lo, (s3 + 4) * delta * scale + lo))

  # --- Spine round application (modular) ---
  def apply_spine_round(current_T, starts, do_interleave=False, reverse_block_parity=False, K=1, add_cross4=False):
    lo, hi, delta = span_delta(current_T)
    blocks = []
    for i, s in enumerate(starts):
      base = s * delta * K - lo
      src = list(reversed(current_T)) if (reverse_block_parity and (i % 2 == 1)) else current_T
      block = [(int(l + base), int(r + base)) for (l, r) in src]
      blocks.append(block)
    # ordering policy
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
    append_classic_connectors(S, starts, delta, lo=0, add_cross4=add_cross4, scale=K)
    return S

  # --- Thin sampling for micro-phase ---
  def thin_seed(current_T, max_seed):
    n = len(current_T)
    if n == 0 or max_seed <= 0:
      return []
    step = max(1, n // max_seed)
    return current_T[::step][:max_seed]

  # --- Micro-phase builder (delta2 primary) ---
  def build_micro_round(current_T, round_id, budget, seed, alt_windows=False, scale_K=1, add_cross4=False):
    if budget <= 0 or not current_T:
      return []
    glo = min(l for l, r in current_T)
    ghi = max(r for l, r in current_T)
    G = max(1, ghi - glo)
    # delta2: smaller subscale — we use half-scale, optionally with K multiplier
    delta2 = max(1, (G // 2))
    # Thin seed size tuned relative to budget
    per_block_target = max(8, min(96, max(8, budget // 12)))
    U = thin_seed(current_T, per_block_target)
    if not U:
      return []
    ulo = min(l for l, r in U)
    # Choose deterministic window offsets based on seed and whether alternate windows requested
    s = seed
    # derive small fractional shifts deterministically (0..9)
    shift_idx = (s ^ (round_id * 1237)) % 10
    small_shift = (shift_idx - 5) * 0.01  # -0.05 .. 0.04 in steps of 0.01
    if not alt_windows:
      window_fracs = [
        (0.10 + small_shift, 0.20 + small_shift),
        (0.32 + small_shift, 0.42 + small_shift),
        (0.55 + small_shift, 0.65 + small_shift),
        (0.78 + small_shift, 0.88 + small_shift),
      ]
    else:
      # alternate window set intentionally different
      window_fracs = [
        (0.05, 0.15),
        (0.28, 0.38),
        (0.60, 0.70),
        (0.82, 0.92),
      ]
    # Clamp fracs to [0.03,0.97] and ensure ordering
    wf = []
    for a, b in window_fracs:
      a = max(0.03, min(0.97, a))
      b = max(a + 0.01, min(0.98, b))
      wf.append((a, b))
    # Build micro-blocks aligned to windows
    starts = template_bank[0]  # keep the classic starts for micro blocks
    blocks = []
    for (fa, fb) in wf:
      win_lo = glo + int(round(fa * G))
      base = win_lo - ulo
      block = [(int(l + base), int(r + base)) for (l, r) in U]
      # break symmetry deterministically
      tag = (seed + int(round(fa * 100))) % 3
      if tag % 2 == 1:
        block = list(reversed(block))
      blocks.append(block)
    # Interleave blocks with parity determined by round_id
    micro = []
    maxlen = max(len(b) for b in blocks)
    order = list(range(len(blocks)))
    if (round_id + seed) % 2 == 1:
      order.reverse()
    for i in range(maxlen):
      for idx in order:
        blk = blocks[idx]
        if i < len(blk):
          micro.append(blk[i])
    # Add localized connectors at delta2 scale (add_cross4 optionally)
    append_classic_connectors(micro, starts, delta2, lo=glo, add_cross4=add_cross4, scale=1)
    # Add a few sparse caps to press FF locally
    capA = (glo + max(1, delta2 // 5), glo + max(1, int(1.6 * delta2)))
    capB = (glo + max(1, int(0.9 * delta2)), glo + max(2, int(2.4 * delta2)))
    mid = glo + G // 2
    capC = (mid - max(1, delta2 // 10), mid + max(1, delta2 // 10))
    for cap in (capA, capB, capC):
      if cap[1] > cap[0]:
        micro.append(cap)
    # Trim micro to budget
    if len(micro) > budget:
      micro = micro[:budget]
    return micro

  # --- MAIN construction ---
  # Seed with one unit interval for classic KT growth
  T = [(0, 1)]
  # Determine safe number of spine rounds given CAP
  rounds, _ = max_rounds_within_cap(len(T), max_spine_rounds)
  # Use at most max_spine_rounds but stop earlier if cap risk
  rounds = min(rounds, max_spine_rounds)
  # apply spine rounds (rotating templates, parity interleaving rules, K density)
  seed = base_seed
  for ridx in range(rounds):
    starts = template_bank[ridx % len(template_bank)]
    # parity choices to diversify order
    do_interleave = (ridx % 2 == 0)  # interleave on even rounds
    reverse_block_parity = (ridx % 2 == 1)
    add_cross4 = (ridx == rounds - 1)  # only final spine round gets a cross4 to avoid omega blow-up
    # apply denser packing for mid-late rounds optionally
    K = use_density_K if ridx >= 1 else 1
    T = apply_spine_round(T, starts, do_interleave=do_interleave, reverse_block_parity=reverse_block_parity, K=K, add_cross4=add_cross4)
    seed = lcg_next(seed)

    # Capacity guard: if we are dangerously close, stop adding spine rounds
    if len(T) >= CAP - 64:
      break

  # Early return if near capacity
  if len(T) >= CAP - 16:
    return T[:CAP]

  # --- Micro-phase A: primary delta2 micro-round (budgeted) ---
  remaining = CAP - len(T)
  if remaining > 16:
    # derive a seed for micro-phase (deterministic)
    seedA = lcg_next(seed)
    # allocate a modest share
    budgetA = max(8, remaining // 3)
    microA = build_micro_round(T, round_id=0, budget=budgetA, seed=seedA, alt_windows=False, scale_K=1, add_cross4=True)
    if microA:
      room = CAP - len(T)
      if len(microA) > room:
        microA = microA[:room]
      T.extend(microA)
    seed = lcg_next(seedA)

  # --- Micro-phase B: alternate delta2 micro-round with different windows and seed ---
  remaining = CAP - len(T)
  if remaining > 12:
    seedB = lcg_next(seed)
    budgetB = max(8, remaining // 2)
    microB = build_micro_round(T, round_id=1, budget=budgetB, seed=seedB, alt_windows=True, scale_K=1, add_cross4=False)
    if microB:
      room = CAP - len(T)
      if len(microB) > room:
        microB = microB[:room]
      T.extend(microB)
    seed = lcg_next(seedB)

  # --- Final micro-tail to hit subtle late FF interactions (tiny) ---
  remaining = CAP - len(T)
  if remaining > 0:
    # Insert a few long caps near the tail to intersect many active colors
    lo, hi, delta = span_delta(T)
    d2 = max(1, delta // 3)
    tail_caps = [
      (lo + 1 * d2, lo + 6 * d2),
      (lo + 3 * d2, lo + 9 * d2),
      (hi - 7 * d2, hi - 2 * d2),
    ]
    # insert staggered near the end to maximize intersection with many active colors
    for i, iv in enumerate(tail_caps):
      if remaining <= 0:
        break
      pos = len(T) - (i * 2 + 1)
      if pos < 0:
        T.append(iv)
      else:
        T.insert(pos, iv)
      remaining -= 1

  # enforce final cap
  if len(T) > CAP:
    T = T[:CAP]
  return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()