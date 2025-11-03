# EVOLVE-BLOCK-START

def construct_intervals(seed_count=1):
  """
  Construct a sequence of intervals on the real line, in the order presented to FirstFit,
  to maximize FF colors divided by the clique number.
  Returns:
    intervals: list of (l, r) tuples (open intervals).
  """

  # Global capacity and guards
  CAP = 9800
  CAP_MARGIN = 32         # guard margin for micro-phases
  BASE_SEED = 123456      # deterministic seed chain
  SPINE_ROUNDS = 6
  SPINE_DENSITY_K = 1     # exploratory density multiplier (kept at 1 by default)
  LR_CONNECTORS = True    # long-range connectors, CAP-gated and conservative
  PHASE2_ITERS = 2        # request up to two micro iterations (A then B)

  # Rotating four-start templates for the KT spine
  template_bank = [
    (2, 6, 10, 14),  # classic KT
    (1, 5, 9, 13),   # left-shifted
    (3, 7, 11, 15),  # right-shifted
    (4, 8, 12, 16),  # stretched-right
  ]

  # Seed with a single unit interval; multi-seed tends to inflate omega too early
  T = [(0, 1)] if seed_count <= 1 else [(0, 1)]

  # ------------------------ helpers ------------------------
  def _span_delta(current_T):
    lo = min(l for l, r in current_T)
    hi = max(r for l, r in current_T)
    delta = hi - lo
    if delta <= 0:
      delta = 1
    return lo, hi, delta

  def _append_connectors(S, starts, delta):
    s0, s1, s2, s3 = starts
    S.append(((s0 - 1) * delta, (s1 - 1) * delta))  # left cap
    S.append(((s2 + 2) * delta, (s3 + 2) * delta))  # right cap
    S.append(((s0 + 2) * delta, (s2 - 1) * delta))  # cross 1
    S.append(((s1 + 2) * delta, (s3 - 1) * delta))  # cross 2

  def _apply_round(current_T, starts, do_interleave=False, reverse_order=False):
    lo, hi, delta = _span_delta(current_T)
    # Build translated blocks
    blocks = []
    for s in starts:
      base = s * delta - lo
      block = [(l + base, r + base) for (l, r) in current_T]
      blocks.append(block)

    # Assemble with optional interleaving and reverse order
    S = []
    if do_interleave:
      maxlen = max(len(b) for b in blocks)
      order = list(range(len(blocks)))
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

    # Classic connectors
    _append_connectors(S, starts, delta)
    return S

  def _predict_next_size(sz):
    # KT growth per round: size -> 4*size + 4 (conservative bound)
    return 4 * sz + 4

  def _det_shift(seed, idx, scale=0.01, mod=7):
    # Tiny, deterministic fractional shift in [-scale, +scale]
    h = (seed ^ (idx * 0x9E3779B1)) & 0xFFFFFFFF
    t = (h % (2 * mod + 1)) - mod
    return (t / float(mod)) * scale

  def _guard_append(T, items, cap):
    if not items:
      return T
    room = cap - len(T)
    if room <= 0:
      return T
    if len(items) > room:
      items = items[:room]
    T.extend(items)
    return T

  # ------------------------ Stage 1: KT spine ------------------------
  # Alternating interleaving and reverse block order; rotate templates each round
  for ridx in range(SPINE_ROUNDS):
    # Guard next round against CAP, reserving space for micro-phases (CAP_MARGIN)
    if _predict_next_size(len(T)) > CAP - CAP_MARGIN:
      break
    starts = template_bank[ridx % len(template_bank)]
    do_inter = (ridx % 2 == 0)
    rev = (ridx % 2 == 1)
    T = _apply_round(T, starts, do_interleave=do_inter, reverse_order=rev)

    # Optional density multiplier K: add a second pass with a different template
    # CAP-guarded; by default disabled (K=1) to preserve proven strong baseline.
    if SPINE_DENSITY_K >= 2 and len(T) < CAP - CAP_MARGIN:
      # Choose a complementary start template deterministically
      starts2 = template_bank[(ridx + 1) % len(template_bank)]
      # Predict size of adding another round’s worth; bail if it risks CAP
      if _predict_next_size(len(T)) <= CAP - CAP_MARGIN:
        T = _apply_round(T, starts2, do_interleave=not do_inter, reverse_order=not rev)

  # Conservative long-range connectors (cross-scale), CAP-gated
  if LR_CONNECTORS and len(T) <= CAP - CAP_MARGIN // 2:
    lo, hi, delta = _span_delta(T)
    # Link distant starts; use the canonical (2, 6, 10, 14) anchors
    s0, s1, s2, s3 = (2, 6, 10, 14)
    lr = [
      (lo + (s0 + 4) * delta, lo + (s3 + 4) * delta),
      (lo + (s1 + 5) * delta, lo + (s2 + 5) * delta),
    ]
    lr = [(a, b) for (a, b) in lr if b > a]
    T = _guard_append(T, lr, CAP)

  # If nearly at capacity, return
  if len(T) >= CAP - CAP_MARGIN:
    return T

  # ------------------------ Stage 2: dual micro-phases ------------------------
  # Micro-phase builder: thin evenly-spaced seed lifted into fractional windows.
  def build_micro(current_T, budget, iter_id=0, alt=False):
    if not current_T or budget <= 8:
      return []

    glo = min(l for l, r in current_T)
    ghi = max(r for l, r in current_T)
    G = max(1, ghi - glo)

    # Thin, evenly spaced seed (bounded)
    seed_sz = max(8, min(40 if alt else 32, len(current_T) // (250 if alt else 300)))
    stride = max(1, len(current_T) // max(1, seed_sz))
    U = [current_T[i] for i in range(0, len(current_T), stride)][:seed_sz]
    if not U:
      return []

    ulo = min(l for l, r in U)

    # Window families with tiny deterministic shifts derived from BASE_SEED
    if not alt:
      # Family A
      window_fracs = [(0.12, 0.22), (0.35, 0.45), (0.58, 0.68), (0.80, 0.90)]
    else:
      # Family B
      window_fracs = [(0.05, 0.15), (0.28, 0.38), (0.60, 0.70), (0.82, 0.92)]

    # Apply tiny deterministic shifts, clipped to [0.03, 0.95]
    fracs = []
    for k, (fa, fb) in enumerate(window_fracs):
      sh = _det_shift(BASE_SEED + (17 if alt else 0), iter_id * 7 + k, scale=0.006, mod=5)
      a = max(0.03, min(0.92, fa + sh))
      b = max(a + 0.01, min(0.95, fb + sh))
      fracs.append((a, b))

    # Build translated micro-blocks aligned to windows
    blocks = []
    for k, (fa, fb) in enumerate(fracs):
      win_lo = glo + int(round(fa * G))
      base = win_lo - ulo
      block = [(l + base, r + base) for (l, r) in U]
      # Deterministic reversal pattern to break symmetry
      if ((k + iter_id) % 2 == 1) ^ alt:
        block = list(reversed(block))
      blocks.append(block)

    # Interleave blocks; alternate order between families
    micro = []
    maxlen = max(len(b) for b in blocks)
    order = list(range(len(blocks)))
    if alt:
      order.reverse()
    for i in range(maxlen):
      for idx in order:
        blk = blocks[idx]
        if i < len(blk):
          micro.append(blk[i])

    # Fractional-span connectors
    connectors = [
      (glo + int(round(0.08 * G)), glo + int(round(0.30 * G))),  # left cap
      (glo + int(round(0.60 * G)), glo + int(round(0.92 * G))),  # right cap
      (glo + int(round(0.26 * G)), glo + int(round(0.56 * G))),  # cross 1
      (glo + int(round(0.44 * G)), glo + int(round(0.78 * G))),  # cross 2
    ]
    if alt:
      # One long cross in the alternate family
      connectors.append((glo + int(round(0.18 * G)), glo + int(round(0.84 * G))))
    micro.extend([(a, b) for (a, b) in connectors if b > a])

    # CAP trim
    if len(micro) > budget:
      micro = micro[:budget]
    return micro

  # Calibrated execution of micro-phases using CAP guard
  # Execute at most PHASE2_ITERS micro-steps: first A, then B
  # Each step requires len(T) <= CAP - CAP_MARGIN
  steps = max(0, int(PHASE2_ITERS))
  if steps >= 1 and len(T) <= CAP - CAP_MARGIN:
    room = CAP - len(T)
    microA = build_micro(T, room, iter_id=0, alt=False)
    T = _guard_append(T, microA, CAP)

  if steps >= 2 and len(T) <= CAP - CAP_MARGIN:
    room = CAP - len(T)
    microB = build_micro(T, room, iter_id=1, alt=True)
    T = _guard_append(T, microB, CAP)

  # A tiny CAP-safe long-range layer after micro-phase (optional)
  if LR_CONNECTORS and len(T) <= CAP - 2:
    lo, hi, delta = _span_delta(T)
    s0, s3 = 2, 14
    tail_lr = [
      (lo + (s0 + 6) * delta, lo + (s3 + 6) * delta),
      (lo + (s0 + 5) * delta, lo + (s3 + 7) * delta),
    ]
    tail_lr = [(a, b) for (a, b) in tail_lr if b > a]
    T = _guard_append(T, tail_lr, CAP)

  # Final capacity trim (safety)
  if len(T) > CAP:
    T = T[:CAP]
  return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()