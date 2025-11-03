# EVOLVE-BLOCK-START

def construct_intervals(rounds=6,
                        rotate_starts=True,
                        reverse_block_parity=True,
                        interleave_blocks=True,
                        phase2_iters=1,
                        cross4_enabled=True):
  """
  Deterministic KT-style rotating spine with CAP-aware densifier,
  parity interleaving, two micro-phases, and long connectors.

  Returns:
    intervals: list[(l, r)] of open intervals in FirstFit presentation order.
  """

  # Hard capacity guard to keep total intervals < 10000
  CAP = 9800

  # Deterministic seed chain to unify all choices
  BASE_SEED = 0x9E3779B97F4A7C15

  # Spine start templates
  spine_starts = (2, 6, 10, 14)
  template_bank = [
    (2, 6, 10, 14),  # classic KT
    (1, 5, 9, 13),   # left-shifted
    (3, 7, 11, 15),  # right-shifted
    (4, 8, 12, 16),  # stretched-right
  ]

  # CAP-aware spine density multiplier ~1.25 (adds one auxiliary block on some rounds)
  SPINE_DENSITY_K = 1.25

  # Seed
  T = [(0, 1)]

  # Helpers
  def _span_delta(current_T):
    lo = min(l for l, r in current_T)
    hi = max(r for l, r in current_T)
    delta = hi - lo
    if delta <= 0:
      delta = 1
    return lo, hi, delta

  def _append_classic_connectors(S, starts, delta):
    s0, s1, s2, s3 = starts
    S.append(((s0 - 1) * delta, (s1 - 1) * delta))  # left cap
    S.append(((s2 + 2) * delta, (s3 + 2) * delta))  # right cap
    S.append(((s0 + 2) * delta, (s2 - 1) * delta))  # cross 1
    S.append(((s1 + 2) * delta, (s3 - 1) * delta))  # cross 2

  def _apply_round(current_T, starts, do_interleave=False, reverse_order=False,
                   add_cross4=False, enable_densify=False, cap_left=1<<60):
    lo, hi, delta = _span_delta(current_T)

    # Build translated blocks from current_T at the 4 starts
    blocks = []
    for s in starts:
      base = s * delta - lo
      block = [(l + base, r + base) for (l, r) in current_T]
      blocks.append(block)

    # Optional densifier: one mid-block between s1 and s2 on selected rounds (CAP-aware)
    densify_block = []
    if enable_densify:
      s0, s1, s2, s3 = starts
      s_mid = (s1 + s2) // 2  # e.g., 8 for (2,6,10,14)
      base_mid = s_mid * delta - lo
      densify_block = [(l + base_mid, r + base_mid) for (l, r) in current_T]

    # Assemble S with optional interleaving and reverse order
    S = []
    seq_blocks = list(blocks)
    if densify_block and len(densify_block) <= cap_left:
      # Insert densify block between the 2nd and 3rd blocks for better mixing
      seq_blocks = [blocks[0], blocks[1], densify_block, blocks[2], blocks[3]]

    if do_interleave:
      maxlen = max(len(b) for b in seq_blocks)
      order = list(range(len(seq_blocks)))
      if reverse_order:
        order.reverse()
      for i in range(maxlen):
        for idx in order:
          blk = seq_blocks[idx]
          if i < len(blk):
            S.append(blk[i])
    else:
      if reverse_order:
        seq_blocks = list(reversed(seq_blocks))
      for blk in seq_blocks:
        S.extend(blk)

    # Classic connectors
    _append_classic_connectors(S, starts, delta)

    # Optional long-range cross connector inside a single round
    if add_cross4:
      s0, s1, s2, s3 = starts
      L = (s0 + 4) * delta
      R = (s3 + 4) * delta
      if R > L:
        S.append((L, R))

    return S, delta

  # Backbone: rotate templates; interleave on even rounds, sequential on odd;
  # reverse block order on odd if reverse_block_parity is True
  last_starts = spine_starts
  last_delta = 1
  for ridx in range(max(0, int(rounds))):
    # Predict size for safety: sz -> ~4*sz + 4 connectors (+ optional densify block)
    connectors_count = 4 + (1 if (cross4_enabled and ridx >= 4) else 0)
    projected = 4 * len(T) + connectors_count
    densify_ok = (SPINE_DENSITY_K > 1.0) and (ridx % 2 == 0) and (len(T) * 5 + connectors_count <= (CAP - len(T)))
    if densify_ok:
      projected += len(T)
    if projected > CAP:
      break

    # Decide starts
    starts = template_bank[ridx % len(template_bank)] if rotate_starts else spine_starts
    last_starts = starts

    # Parity schedule: interleave on even rounds only
    do_inter = bool(interleave_blocks and (ridx % 2 == 0))
    rev = bool(reverse_block_parity and (ridx % 2 == 1))

    # Deterministic toggle for densifier per round via seed
    per_round_seed = (BASE_SEED ^ (ridx * 0xBF58476D1CE4E5B9)) & ((1 << 61) - 1)
    densify_flag = densify_ok and ((per_round_seed >> 7) % 3 != 0)  # activate 2/3 of eligible even rounds

    # Long cross inside round: only on late rounds to keep omega modest
    add_c4 = bool(cross4_enabled and ridx >= 4 and (per_round_seed & 1) == 0)

    cap_left = CAP - len(T)
    S, used_delta = _apply_round(
      T, starts,
      do_interleave=do_inter,
      reverse_order=rev,
      add_cross4=add_c4,
      enable_densify=densify_flag,
      cap_left=cap_left
    )
    last_delta = used_delta

    if len(S) > CAP:
      S = S[:CAP]
    T = S

    if len(T) >= CAP:
      T = T[:CAP]
      return T

  # Early return if near capacity
  if len(T) >= CAP - 8:
    return T[:CAP]

  # Micro-phase A: small set of long caps near the end to boost mixing
  lo, hi, G = _span_delta(T)
  def cap_at(a_frac, b_frac):
    L = lo + max(1, int(round(a_frac * G)))
    R = lo + max(1, int(round(b_frac * G)))
    if R <= L:
      R = L + 1
    return (L, R)
  tail_caps = [cap_at(0.08, 0.60), cap_at(0.25, 0.75), cap_at(0.75, 0.92)]
  for c in tail_caps:
    if len(T) >= CAP:
      break
    T.append(c)

  if len(T) >= CAP - 16:
    return T[:CAP]

  # Micro-phase builder: two deterministic window families, thin seed, interleaving
  def build_micro_round(current_T, budget, family_id=0):
    if not current_T or budget <= 8:
      return []

    glo = min(l for l, r in current_T)
    ghi = max(r for l, r in current_T)
    span = max(1, ghi - glo)

    # Thin evenly-spaced seed
    n = len(current_T)
    seed_sz = max(8, min(36, n // 280))
    stride = max(1, n // max(1, seed_sz))
    # deterministic offset for diversity by family
    off = ((BASE_SEED ^ (family_id * 0x94D049BB133111EB)) % max(1, stride))
    U = [current_T[(i + off) % n] for i in range(0, n, stride)][:seed_sz]
    if not U:
      return []
    ulo = min(l for l, r in U)

    # Two window families
    if family_id == 0:
      windows = [(0.12, 0.22), (0.35, 0.45), (0.58, 0.68), (0.80, 0.90)]
    else:
      windows = [(0.06, 0.14), (0.28, 0.38), (0.54, 0.64), (0.78, 0.88)]

    # Build translated micro-blocks
    blocks = []
    for (fa, fb) in windows:
      fa = max(0.05, min(0.90, fa))
      fb = max(fa + 0.02, min(0.95, fb))
      win_lo = glo + int(round(fa * span))
      base = win_lo - ulo
      block = [(l + base, r + base) for (l, r) in U]
      # Alternate internal reversal deterministically
      if int((fa * 100)) % 2 == (family_id % 2):
        block = list(reversed(block))
      blocks.append(block)

    # Interleave blocks; reverse order for family 1
    micro = []
    order = list(range(len(blocks)))
    if family_id % 2 == 1:
      order.reverse()
    maxlen = max(len(b) for b in blocks)
    for i in range(maxlen):
      for idx in order:
        blk = blocks[idx]
        if i < len(blk):
          micro.append(blk[i])

    # Fractional-span micro connectors
    connectors = [
      (glo + int(round(0.08 * span)), glo + int(round(0.30 * span))),
      (glo + int(round(0.60 * span)), glo + int(round(0.92 * span))),
      (glo + int(round(0.26 * span)), glo + int(round(0.56 * span))),
      (glo + int(round(0.44 * span)), glo + int(round(0.78 * span))),
    ]
    for a, b in connectors:
      if b > a:
        micro.append((a, b))

    if len(micro) > budget:
      micro = micro[:budget]
    return micro

  # Execute micro-phase(s) with capacity guard
  steps = min(max(0, int(phase2_iters)), 2)
  if steps >= 1:
    room = CAP - len(T)
    if room > 8:
      microA = build_micro_round(T, room, family_id=0)
      if microA:
        add = microA[:room]
        T.extend(add)

  if steps >= 1:
    room = CAP - len(T)
    if room > 8:
      microB = build_micro_round(T, room, family_id=1)
      if microB:
        add = microB[:room]
        T.extend(add)

  # Long-range deterministic cross-scale connectors (4–6 delta) anchored on last spine geometry
  if cross4_enabled:
    room = CAP - len(T)
    if room > 0:
      d = max(1, int(last_delta))
      s0, s1, s2, s3 = last_starts
      long_pairs = [
        ((s0 - 2) * d, (s2 + 3) * d),
        ((s1 - 1) * d, (s3 + 4) * d),
      ]
      for a, b in long_pairs:
        if len(T) >= CAP:
          break
        if b > a:
          T.append((a, b))

  # Final capacity trim
  if len(T) > CAP:
    T = T[:CAP]
  return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()