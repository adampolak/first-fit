# EVOLVE-BLOCK-START

def construct_intervals(rounds=6,
                        rotate_starts=True,
                        reverse_block_parity=True,
                        interleave_blocks=True,
                        phase2_iters=1):
  """
  KT-style spine with early K-densify and a deterministic, two-family micro-phase.

  Parameters:
    rounds (int): requested backbone rounds (cap-aware actual depth).
    rotate_starts (bool): rotate among four start templates for the backbone.
    reverse_block_parity (bool): flip block order on odd rounds in sequential mode.
    interleave_blocks (bool): enable interleaving on even rounds (parity-based).
    phase2_iters (int): retained for compatibility; micro-phase is CAP-driven.

  Returns:
    intervals: list of (l, r) integer tuples (open intervals), in FF arrival order.
  """

  # ----------------------
  # Global capacity policy
  # ----------------------
  CAP = 9800
  MICRO_RESERVE = 6000  # reserve for late micro phases to push FirstFit up

  # ----------------------
  # Backbone templates
  # ----------------------
  spine_starts_default = (2, 6, 10, 14)
  template_bank = [
    (2, 6, 10, 14),  # classic KT
    (1, 5, 9, 13),   # left-shifted
    (3, 7, 11, 15),  # right-shifted
    (4, 8, 12, 16),  # stretched-right
  ]

  # Deterministic densify schedule: add one extra intermediate block on early rounds
  DENSIFY_ROUNDS = {0, 1}           # apply on the first two backbone rounds
  DENSIFY_COUNT_PER_ROUND = 1       # exactly one extra block at midpoint between s1 and s2

  # ----------------------
  # Helpers
  # ----------------------
  def _span_delta(current_T):
    lo = min(l for l, r in current_T)
    hi = max(r for l, r in current_T)
    delta = hi - lo
    if delta <= 0:
      delta = 1
    return lo, hi, int(delta)

  def _append_connectors(S, starts, delta, base_shift=0):
    # Classic four connectors (scaled by delta). base_shift allows localized micro connectors.
    s0, s1, s2, s3 = starts
    S.append((base_shift + (s0 - 1) * delta, base_shift + (s1 - 1) * delta))  # left cap
    S.append((base_shift + (s2 + 2) * delta, base_shift + (s3 + 2) * delta))  # right cap
    S.append((base_shift + (s0 + 2) * delta, base_shift + (s2 - 1) * delta))  # cross 1
    S.append((base_shift + (s1 + 2) * delta, base_shift + (s3 - 1) * delta))  # cross 2

  def _apply_backbone_round(current_T, starts, ridx, do_interleave):
    lo, hi, delta = _span_delta(current_T)

    # Build the four translated blocks
    blocks = []
    for s in starts:
      base = s * delta - lo
      block = [(l + base, r + base) for (l, r) in current_T]
      blocks.append(block)

    # Optional densify: add exactly one extra block at midpoint between starts[1] and starts[2]
    densify = (ridx in DENSIFY_ROUNDS)
    densify_block = None
    if densify and DENSIFY_COUNT_PER_ROUND >= 1:
      mid = starts[1] + (starts[2] - starts[1]) // 2  # integer midpoint
      base_mid = mid * delta - lo
      densify_block = [(l + base_mid, r + base_mid) for (l, r) in current_T]

    # Compose S with interleaving on even rounds; sequential (with optional reverse) on odd
    S = []
    if do_interleave:
      # Interleave four (or five) blocks
      interleave_list = list(blocks)
      if densify_block is not None:
        # place densify block in between second and third to maximize mixing
        interleave_list = [blocks[0], blocks[1], densify_block, blocks[2], blocks[3]]
      maxlen = max(len(b) for b in interleave_list)
      for i in range(maxlen):
        for blk in interleave_list:
          if i < len(blk):
            S.append(blk[i])
    else:
      # Sequential; optionally reverse block order on odd rounds
      seq_blocks = list(blocks)
      if reverse_block_parity and (ridx % 2 == 1):
        seq_blocks = list(reversed(seq_blocks))
      # Insert densify block after the first two blocks to emulate midpoint placement
      if densify_block is not None:
        merged = []
        merged.extend(seq_blocks[0])
        merged.extend(seq_blocks[1])
        merged.extend(densify_block)
        merged.extend(seq_blocks[2])
        merged.extend(seq_blocks[3])
        S.extend(merged)
      else:
        for blk in seq_blocks:
          S.extend(blk)

    # Append the classic connectors at the backbone scale
    _append_connectors(S, starts, delta, base_shift=0)

    return S

  # ----------------------
  # Micro-phase machinery
  # ----------------------
  # Two deterministic families of windows to be applied alternately after the backbone
  window_fracs1 = [(0.12, 0.22), (0.35, 0.45), (0.58, 0.68), (0.80, 0.90)]
  window_fracs2 = [(0.06, 0.14), (0.28, 0.38), (0.54, 0.64), (0.78, 0.88)]

  # Seed sizing parameters (thin seed to control omega)
  SEED_BASE = 32
  SEED_MAX = 64
  SEED_DIV = 280  # seed size roughly |T| // SEED_DIV, clamped in [8, SEED_MAX]

  def _thin_seed(current_T, seed_cap):
    """Evenly spaced thin seed from current_T of size <= seed_cap."""
    n = len(current_T)
    if n == 0 or seed_cap <= 0:
      return []
    step = max(1, n // seed_cap)
    U = current_T[::step][:seed_cap]
    return U

  def _build_micro_blocks_from_windows(current_T, window_fracs, reverse_order=False):
    if not current_T:
      return []

    glo = min(l for l, r in current_T)
    ghi = max(r for l, r in current_T)
    G = max(1, ghi - glo)

    # Seed size adapts to current_T but stays thin
    seed_sz = max(8, min(SEED_MAX, max(SEED_BASE, len(current_T) // SEED_DIV)))
    U = _thin_seed(current_T, seed_sz)
    if not U:
      return []

    ulo = min(l for l, r in U)

    # Build translated micro-blocks
    blocks = []
    for (fa, fb) in window_fracs:
      fa = max(0.05, min(0.90, fa))
      fb = max(0.10, min(0.95, fb))
      win_lo = glo + int(round(fa * G))
      base = win_lo - ulo
      block = [(l + base, r + base) for (l, r) in U]
      # small internal reversal to break symmetry
      if int(round(100 * fa)) // 5 % 2 == 1:
        block = list(reversed(block))
      blocks.append(block)

    # Interleave blocks; optionally reverse block order
    micro = []
    order = list(range(len(blocks)))
    if reverse_order:
      order.reverse()
    maxlen = max(len(b) for b in blocks)
    for i in range(maxlen):
      for idx in order:
        blk = blocks[idx]
        if i < len(blk):
          micro.append(blk[i])

    # Fractional-span connectors at micro scale (4–6 delta-span), anchored to glo
    micro_connectors = [
      (glo + int(round(0.08 * G)), glo + int(round(0.30 * G))),  # left cap
      (glo + int(round(0.60 * G)), glo + int(round(0.92 * G))),  # right cap
      (glo + int(round(0.26 * G)), glo + int(round(0.56 * G))),  # cross 1
      (glo + int(round(0.44 * G)), glo + int(round(0.78 * G))),  # cross 2
      (glo + int(round(0.12 * G)), glo + int(round(0.64 * G))),  # longer cross
      (glo + int(round(0.36 * G)), glo + int(round(0.88 * G))),  # longer cross
    ]
    for a, b in micro_connectors:
      if b > a:
        micro.append((a, b))

    return micro

  # ----------------------
  # Backbone planning under CAP with densify rounds early
  # ----------------------
  # Simulate size evolution to determine allowable rounds with micro reserve
  def _simulate_backbone_size(base_size, req_rounds):
    size = base_size
    rounds_run = 0
    for ridx in range(max(0, int(req_rounds))):
      copies = 4 + (1 if ridx in DENSIFY_ROUNDS else 0)
      next_size = copies * size + 4  # +4 connectors
      # Stop if exceeding cap - reserve for micro phase
      if next_size > CAP - MICRO_RESERVE:
        break
      size = next_size
      rounds_run += 1
    return rounds_run, size

  # ----------------------
  # Construct the intervals
  # ----------------------
  T = [(0, 1)]  # single seed to keep omega low early

  # Stage 1: backbone with early densify and parity-based interleaving
  depth, _ = _simulate_backbone_size(len(T), rounds)
  for ridx in range(depth):
    starts = template_bank[ridx % len(template_bank)] if rotate_starts else spine_starts_default
    do_inter = bool(interleave_blocks and (ridx % 2 == 0))
    T = _apply_backbone_round(T, starts, ridx, do_inter)
    if len(T) >= CAP - MICRO_RESERVE // 2:
      # Safety break in pathological cases
      break

  # If the backbone already consumed most of the budget, return it
  if len(T) >= CAP - 32:
    return T[:CAP]

  # Stage 2: deterministic two-family micro-phase loop until CAP
  # Alternate between window_fracs1 (forward interleave) and window_fracs2 (reverse interleave)
  family_toggle = 0
  while len(T) < CAP:
    remaining = CAP - len(T)
    if remaining <= 0:
      break

    if family_toggle % 2 == 0:
      micro = _build_micro_blocks_from_windows(T, window_fracs1, reverse_order=False)
    else:
      micro = _build_micro_blocks_from_windows(T, window_fracs2, reverse_order=True)

    if not micro:
      break

    # Trim to remaining capacity with a modest per-iteration limit to keep steps balanced
    step_limit = max(128, min(1200, remaining))
    if len(micro) > step_limit:
      micro = micro[:step_limit]

    T.extend(micro)
    family_toggle += 1

  # Final CAP trim
  if len(T) > CAP:
    T = T[:CAP]
  return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()