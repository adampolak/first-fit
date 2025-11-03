# EVOLVE-BLOCK-START

def construct_intervals(seed_count=1):
  """
  Construct a sequence of intervals on the real line (open intervals),
  presented in FirstFit order, to maximize FF colors while keeping omega small.

  Inputs/Outputs:
    seed_count (int): retains compatibility; builds the same return type
    Returns: list of (l, r) integer tuples representing open intervals (l, r)
  """

  # Global capacity and safety margin
  CAP = 9800
  CAP_MARGIN = 32  # tunable guard to gate micro-phases

  # Backbone configuration parameters
  ROUNDS = 6
  rotate_starts = True
  reverse_block_parity = True
  interleave_blocks = True

  # Deterministic seed chain for small shifts/parity tweaks
  BASE_SEED = 1337

  # Spine density multiplier (K=2 appends an alternate translated set per round)
  K = 2

  # Start templates (rotated deterministically across rounds)
  template_bank = [
    (2, 6, 10, 14),  # classic KT
    (1, 5, 9, 13),   # left-shifted
    (3, 7, 11, 15),  # right-shifted
    (4, 8, 12, 16),  # stretched-right
  ]

  # Seed: single interval, unless user requests multiple (kept deterministic)
  if seed_count <= 1:
    T = [(0, 1)]
  else:
    step = 3
    T = [(i * step, i * step + 1) for i in range(seed_count)]

  # Utility helpers
  def _span_delta(current_T):
    lo = min(l for l, r in current_T)
    hi = max(r for l, r in current_T)
    delta = hi - lo
    if delta <= 0:
      delta = 1
    return lo, hi, delta

  def _append_connectors(S, starts, delta):
    s0, s1, s2, s3 = starts
    S.append(( (s0 - 1) * delta, (s1 - 1) * delta ))  # left cap
    S.append(( (s2 + 2) * delta, (s3 + 2) * delta ))  # right cap
    S.append(( (s0 + 2) * delta, (s2 - 1) * delta ))  # cross 1
    S.append(( (s1 + 2) * delta, (s3 - 1) * delta ))  # cross 2

  def _append_long_range_connectors(S, starts, delta):
    # Sparse long connectors spanning 4–6 deltas; helps couple distant colors
    s0, s1, s2, s3 = starts
    L1 = ( (s0 - 2) * delta, (s3 + 3) * delta )
    L2 = ( (s0 + 1) * delta, (s2 + 3) * delta )
    if L1[1] > L1[0]:
      S.append(L1)
    if L2[1] > L2[0]:
      S.append(L2)

  def _det_round_seed(base, ridx):
    # Simple LCG to derive deterministic per-round "seed"
    return (base * 48271 + 0x9e3779b97f4a7c15 + ridx * 1664525) & 0xFFFFFFFF

  def _apply_round(current_T, starts, ridx, enable_interleave, enable_reverse, k_density=1):
    lo, hi, delta = _span_delta(current_T)
    # Build primary translated blocks
    blocks = []
    for s in starts:
      base = s * delta - lo
      block = [(l + base, r + base) for (l, r) in current_T]
      blocks.append(block)

    # Optional K=2 density: append alternate template shifted by +1
    alt_blocks = []
    if k_density >= 2:
      alt_starts = tuple(x + 1 for x in starts)
      for s in alt_starts:
        base = s * delta - lo
        block = [(l + base, r + base) for (l, r) in current_T]
        # reverse internal order for diversity
        block = list(reversed(block))
        alt_blocks.append(block)

    # Merge blocks according to interleave/reversal policies
    S = []
    sequence = blocks + alt_blocks if alt_blocks else blocks
    if enable_reverse:
      sequence = list(reversed(sequence))

    if enable_interleave and sequence:
      maxlen = max(len(b) for b in sequence)
      for i in range(maxlen):
        for blk in sequence:
          if i < len(blk):
            S.append(blk[i])
    else:
      for blk in sequence:
        S.extend(blk)

    # Append connectors
    _append_connectors(S, starts, delta)
    if alt_blocks:
      _append_connectors(S, tuple(x + 1 for x in starts), delta)

    # Sparse long-range connectors (cross-scale)
    _append_long_range_connectors(S, starts, delta)
    return S

  # Predictive size estimation to keep inside CAP when using K=2
  def _estimate_next_size(sz, k_density):
    # Rough model: per-round cloning into ~4*k blocks plus ~10 connectors per template
    clones = 4 * k_density
    conn = 10 * k_density
    return clones * sz + conn

  # Stage 1: KT-style spine with rotating templates and K-density (guarded by CAP estimate)
  for ridx in range(ROUNDS):
    starts = template_bank[ridx % len(template_bank)] if rotate_starts else (2, 6, 10, 14)
    # Decide interleave and reversal
    do_inter = bool(interleave_blocks and (ridx % 2 == 0))
    do_rev = bool(reverse_block_parity and (ridx % 2 == 1))

    # CAP-aware density choice
    est_next = _estimate_next_size(len(T), K)
    k_eff = K if est_next <= CAP else 1

    T = _apply_round(T, starts, ridx, enable_interleave=do_inter, enable_reverse=do_rev, k_density=k_eff)
    if len(T) >= CAP - CAP_MARGIN:
      # Stop early if too close to capacity
      break

  if len(T) >= CAP - CAP_MARGIN:
    return T

  # Micro-phase A: three long-range caps inserted near the end
  lo, hi, delta = _span_delta(T)
  def cap_at(a_frac, b_frac):
    L = lo + max(1, int(round(a_frac * delta)))
    R = lo + max(1, int(round(b_frac * delta)))
    if R <= L:
      R = L + 1
    return (L, R)
  tail_caps = [cap_at(0.08, 0.60), cap_at(0.25, 0.75), cap_at(0.75, 0.92)]
  # Insert toward the tail to mix ongoing block colors
  def _insert_near_tail(seq, intervals):
    out = list(seq)
    for i, iv in enumerate(intervals):
      pos = len(out) - (i * 3 + 1)
      if pos < 0:
        out.append(iv)
      else:
        out.insert(pos, iv)
    return out
  room = CAP - len(T)
  if room > 0:
    T = _insert_near_tail(T, tail_caps[:room])

  # CAP gate for micro-phases
  if len(T) > CAP - CAP_MARGIN:
    return T

  # Micro-phase builder (two distinct window sets, deterministic shifts)
  def build_micro_round(current_T, budget, iter_id=0, win_set_id=0, seed=0):
    if not current_T or budget <= 8:
      return []

    glo = min(l for l, r in current_T)
    ghi = max(r for l, r in current_T)
    G = max(1, ghi - glo)

    # Thin seed: evenly spaced subsample (deterministic)
    if win_set_id == 0:
      seed_min, seed_max = 8, 40
    else:
      seed_min, seed_max = 8, 32
    seed_sz = max(seed_min, min(seed_max, len(current_T) // 240))
    stride = max(1, len(current_T) // max(1, seed_sz))
    U = [current_T[i] for i in range(0, len(current_T), stride)][:seed_sz]
    if not U:
      return []

    ulo = min(l for l, r in U)

    # Window sets
    W1 = [(0.12, 0.22), (0.35, 0.45), (0.58, 0.68), (0.80, 0.90)]
    W2 = [(0.05, 0.15), (0.28, 0.38), (0.60, 0.70), (0.82, 0.92)]
    windows = W1 if win_set_id == 0 else W2

    # Deterministic small shift from seed/iter to break symmetry (kept inside [0,0.02])
    shift_units = ((seed ^ (iter_id * 2654435761)) & 7)  # 0..7
    tiny_shift = 0.0025 * shift_units

    # Build translated micro blocks aligned to the windows
    blocks = []
    for (fa, fb) in windows:
      fa2, fb2 = fa + tiny_shift, fb + tiny_shift
      fa2 = max(0.04, min(0.91, fa2))
      fb2 = max(fa2 + 0.03, min(0.96, fb2))
      win_lo = glo + int(round(fa2 * G))
      base = win_lo - ulo
      block = [(l + base, r + base) for (l, r) in U]
      # Alternate internal reversal based on window index and iteration
      idx_tag = int(round(fa * 100))
      if ((idx_tag // 5) + iter_id + win_set_id) % 2 == 1:
        block = list(reversed(block))
      blocks.append(block)

    # Interleave blocks (forward on even iter, reverse on odd)
    micro = []
    maxlen = max(len(b) for b in blocks)
    order = list(range(len(blocks)))
    if iter_id % 2 == 1:
      order.reverse()
    for i in range(maxlen):
      for j in order:
        blk = blocks[j]
        if i < len(blk):
          micro.append(blk[i])

    # Fractional-span connectors
    connectors = [
      (glo + int(round(0.08 * G)), glo + int(round(0.30 * G))),  # left cap
      (glo + int(round(0.60 * G)), glo + int(round(0.92 * G))),  # right cap
      (glo + int(round(0.26 * G)), glo + int(round(0.56 * G))),  # cross 1
      (glo + int(round(0.44 * G)), glo + int(round(0.78 * G))),  # cross 2
    ]
    # Add a single long cross for the second window family
    if win_set_id == 1:
      connectors.append((glo + int(round(0.18 * G)), glo + int(round(0.84 * G))))
    for a, b in connectors:
      if b > a:
        micro.append((a, b))

    # Budget guard
    if len(micro) > budget:
      micro = micro[:budget]
    return micro

  # Execute micro-phase A (up to two iterations), capacity-guarded with CAP_MARGIN
  round_seed = _det_round_seed(BASE_SEED, 0)
  stepsA = 2
  for it in range(stepsA):
    room = CAP - len(T)
    if room <= CAP_MARGIN:
      break
    microA = build_micro_round(T, room, iter_id=it, win_set_id=0, seed=round_seed)
    if not microA:
      break
    if len(microA) > room:
      microA = microA[:room]
    T.extend(microA)

  # Execute micro-phase B (one iteration with alternate windows)
  room = CAP - len(T)
  if room > CAP_MARGIN:
    microB = build_micro_round(T, room, iter_id=0, win_set_id=1, seed=_det_round_seed(BASE_SEED, 1))
    if microB:
      if len(microB) > room:
        microB = microB[:room]
      T.extend(microB)

  return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()