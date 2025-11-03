# EVOLVE-BLOCK-START

def construct_intervals(enable_alt_microphase=True):
  """
  Construct a sequence of intervals presented to FirstFit.

  Enhancements:
    - 6-template deterministic backbone (more mixing variety).
    - Controlled interleaving (2/3 of rounds) and parity reversal.
    - Per-round lightweight densification in middle rounds.
    - Deterministic long-range backbone connectors after the spine.
    - Two-phase micro-phase budgeting (60/40 split) with small window shifts.
    - Final normalization to integer endpoints.

  Args:
    enable_alt_microphase (bool): enable a second, alternate micro-phase.

  Returns:
    intervals: list of tuples (l, r) representing open intervals.
  """

  # Capacity guard to keep total intervals < 10000
  CAP = 9800

  # Expanded deterministic template bank (six useful templates).
  template_bank = [
    (2, 6, 10, 14),  # classic KT
    (1, 5, 9, 13),   # left-shifted
    (3, 7, 11, 15),  # right-shifted
    (4, 8, 12, 16),  # stretched-right
    (2, 4, 8, 12),   # compressed-left pair
    (3, 5, 9, 13),   # gentle left pack
  ]

  # Seed with a single unit interval
  T = [(0, 1)]

  # --- Helpers ---------------------------------------------------------------
  def _span(current_T):
    lo = min(l for l, r in current_T)
    hi = max(r for l, r in current_T)
    delta = hi - lo
    if delta <= 0:
      delta = 1
    return lo, hi, delta

  def _apply_round(current_T, starts, do_interleave=False, reverse_order=False,
                   densify_count=0, densify_eps=0, add_cross4=False):
    """
    Apply a template round: translate blocks defined by 'starts', optionally interleave,
    optionally add a small number of densified (slightly shrunken) copies, and append
    classic connectors. This is capacity-neutral in the sense it keeps block sizes
    predictable; densify_count is small (3-4) in practice.
    """
    lo, hi, delta = _span(current_T)

    # Build translated blocks
    blocks = []
    for s in starts:
      base = s * delta - lo
      block = [(l + base, r + base) for (l, r) in current_T]
      blocks.append(block)

    # Assemble S with optional interleaving and reverse order
    S = []
    if do_interleave:
      order = list(range(len(blocks)))
      if reverse_order:
        order.reverse()
      maxlen = max(len(b) for b in blocks)
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

    # Lightweight densification: pick a few intervals from the current seed and
    # append slightly shrunk copies to increase FirstFit pressure without raising clique much.
    if densify_count > 0 and densify_eps > 0 and current_T:
      n = len(current_T)
      step = max(1, n // max(1, densify_count))
      picks = [current_T[i] for i in range(0, n, step)][:densify_count]
      for (l, r) in picks:
        if r - l > 2 * densify_eps:
          S.append((l + densify_eps, r - densify_eps))

    # Classic connectors (Figure-4 style) scaled to the current delta.
    ss = list(starts)
    # Ensure we can index four entries safely
    while len(ss) < 4:
      ss.append(ss[-1] + 4)
    s0, s1, s2, s3 = ss[0], ss[1], ss[2], ss[3]
    S.append(((s0 - 1) * delta, (s1 - 1) * delta))  # left cap
    S.append(((s2 + 2) * delta, (s3 + 2) * delta))  # right cap
    S.append(((s0 + 2) * delta, (s2 - 1) * delta))  # cross 1
    S.append(((s1 + 2) * delta, (s3 - 1) * delta))  # cross 2

    # Optional conservative long cross to tie distant blocks (used sparingly).
    if add_cross4:
      S.append(((s0 + 4) * delta, (s3 + 4) * delta))

    return S

  def build_micro_delta_round(current_T, budget, alt=False, iter_id=0):
    """
    Build a micro-phase: take a thin seed U of current_T and create four micro-blocks
    placed into fractional windows across the global span. The iter_id introduces
    small deterministic shifts/order changes between repeated micro rounds.
    """
    if not current_T or budget <= 8:
      return []

    glo = min(l for l, r in current_T)
    ghi = max(r for l, r in current_T)
    G = max(1, ghi - glo)

    # Seed sizing (slightly larger in alternate micro-phase to improve mixing)
    if not alt:
      seed_sz = max(8, min(48, len(current_T) // 220))
    else:
      seed_sz = max(8, min(48, len(current_T) // 200))
    stride = max(1, len(current_T) // max(1, seed_sz))
    U = [current_T[i] for i in range(0, len(current_T), stride)][:seed_sz]
    if not U:
      return []

    ulo = min(l for l, r in U)

    # Window families with small deterministic shift controlled by iter_id
    if not alt:
      shift = (iter_id % 4) * 0.015
      window_fracs = [
        (0.12 + shift, 0.22 + shift),
        (0.35 + shift, 0.45 + shift),
        (0.58 + shift, 0.68 + shift),
        (0.80 + shift, 0.90 + shift),
      ]
    else:
      window_fracs = [(0.05, 0.15), (0.28, 0.38), (0.60, 0.70), (0.82, 0.92)]

    # Build translated micro-blocks aligned to windows (and break symmetry slightly).
    blocks = []
    for idx, (fa, fb) in enumerate(window_fracs):
      win_lo = glo + int(round(fa * G))
      base = win_lo - ulo
      block = [(l + base, r + base) for (l, r) in U]
      # Reverse some blocks deterministically to break regularity
      if ((idx + iter_id) % 2) == 1:
        block = list(reversed(block))
      blocks.append(block)

    # Interleave micro-blocks with small ordering variation
    micro = []
    maxlen = max(len(b) for b in blocks)
    order = list(range(len(blocks)))
    if (iter_id % 2) == 1:
      order.reverse()
    for i in range(maxlen):
      for idx in order:
        blk = blocks[idx]
        if i < len(blk):
          micro.append(blk[i])

    # Fractional-span connectors at the micro scale
    micro_connectors = [
      (glo + int(round(0.08 * G)), glo + int(round(0.30 * G))),  # left cap
      (glo + int(round(0.60 * G)), glo + int(round(0.92 * G))),  # right cap
      (glo + int(round(0.26 * G)), glo + int(round(0.56 * G))),  # cross 1
      (glo + int(round(0.44 * G)), glo + int(round(0.78 * G))),  # cross 2
    ]
    for a, b in micro_connectors:
      if b > a:
        micro.append((a, b))

    # A longer cross for the alternate micro-phase to improve coupling
    if alt:
      a = glo + int(round(0.18 * G))
      b = glo + int(round(0.84 * G))
      if b > a:
        micro.append((a, b))

    # Trim to available budget
    if len(micro) > budget:
      micro = micro[:budget]
    return micro

  # --- Stage 1: KT backbone rounds (capacity-guarded) ------------------------
  rounds = 6
  for ridx in range(rounds):
    # Predict next size: sz -> 4*sz + 4
    if 4 * len(T) + 4 > CAP:
      break

    starts = template_bank[ridx % len(template_bank)]

    # Controlled interleaving: do_interleave in 2/3 of rounds to increase mixing.
    do_inter = (ridx % 3 != 1)
    rev = (ridx % 2 == 1)

    # Lightweight densification in middle rounds to raise FF pressure
    lo, hi, delta = _span(T)
    densify_count = 3 if ridx in (2, 3, 4) else 0
    densify_eps = max(1, delta // 1000) if densify_count > 0 else 0

    # Conservative long-cross only on the final backbone round
    add_cross4 = (ridx == rounds - 1)

    T = _apply_round(T, starts, do_interleave=do_inter, reverse_order=rev,
                     densify_count=densify_count, densify_eps=densify_eps,
                     add_cross4=add_cross4)

    # Capacity trim safety
    if len(T) > CAP:
      T = T[:CAP]

  # Early exit if nearly at capacity (before micro-phases)
  if len(T) >= CAP - 8:
    # Normalize endpoints to integers before returning
    min_l = min(l for l, r in T)
    if min_l < 0:
      T = [(l - min_l, r - min_l) for (l, r) in T]
    out = []
    for (l, r) in T[:CAP]:
      li = int(round(l)); ri = int(round(r))
      if ri <= li:
        ri = li + 1
      out.append((li, ri))
    return out

  # --- After backbone: deterministic long-range connectors -------------------
  lo, hi, span = _span(T)
  span = max(1, span)
  backbone_connectors = [
    (lo + max(1, int(round(0.04 * span))), lo + max(2, int(round(0.48 * span)))),
    (lo + max(1, int(round(0.52 * span))), lo + max(2, int(round(0.96 * span)))),
    (lo + max(1, int(round(0.16 * span))), lo + max(2, int(round(0.84 * span)))),
  ]
  for c in backbone_connectors:
    if len(T) >= CAP:
      break
    T.append(c)

  # --- Two-phase micro-phase budgeting anchored to the backbone -------------
  remaining = CAP - len(T)
  if remaining > 8:
    budget1 = max(8, int(round(remaining * 0.60)))
    budget2 = max(0, remaining - budget1)

    # Micro-phase A: forward windows, small iter_id shift
    room = CAP - len(T)
    b1 = min(room, budget1)
    if b1 > 8:
      micro1 = build_micro_delta_round(T, b1, alt=False, iter_id=0)
      if micro1:
        if len(micro1) > b1:
          micro1 = micro1[:b1]
        T.extend(micro1)

    # Recompute room and allocate to alternate micro-phase
    room = CAP - len(T)
    b2 = min(room, budget2)
    if enable_alt_microphase and b2 > 8:
      micro2 = build_micro_delta_round(T, b2, alt=True, iter_id=1)
      if micro2:
        if len(micro2) > b2:
          micro2 = micro2[:b2]
        T.extend(micro2)

  # Final capacity trim and normalization to integer endpoints
  if len(T) > CAP:
    T = T[:CAP]

  if not T:
    return []

  min_l = min(l for l, r in T)
  if min_l < 0:
    T = [(l - min_l, r - min_l) for (l, r) in T]

  intervals = []
  for (l, r) in T[:CAP]:
    li = int(round(l)); ri = int(round(r))
    if ri <= li:
      ri = li + 1
    intervals.append((li, ri))

  return intervals

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()