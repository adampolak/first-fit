# EVOLVE-BLOCK-START

def construct_intervals(rounds=6,
                        rotate_starts=True,
                        reverse_block_parity=True,
                        interleave_blocks=True,
                        phase2_iters=2):
  """
  Rotating four-template KT spine with dual deterministic micro-phases to boost FirstFit.

  Parameters:
    rounds (int): main expansion depth; internally capped for safety and improved ratio.
    rotate_starts (bool): kept for interface compatibility (spine rotation is deterministic).
    reverse_block_parity (bool): kept for interface compatibility.
    interleave_blocks (bool): kept for interface compatibility.
    phase2_iters (int): requested micro iterations; actual use is capacity-guarded (up to two).

  Returns:
    intervals: list of (l, r) integer tuples, open intervals, in FF presentation order.
  """

  CAP = 9800  # hard capacity bound

  # Deterministic start templates (classic and shifted)
  template_bank = [
    (2, 6, 10, 14),  # classic KT
    (1, 5, 9, 13),   # left-shifted
    (3, 7, 11, 15),  # right-shifted
    (4, 8, 12, 16),  # stretched-right
  ]

  # Helpers
  def _span_delta(T):
    lo = min(l for l, _ in T)
    hi = max(r for _, r in T)
    delta = hi - lo
    if delta <= 0:
      delta = 1
    return lo, hi, delta

  def _apply_round(current_T, starts, do_interleave=False, reverse_order=False):
    # Compute span and delta
    lo, hi, delta = _span_delta(current_T)

    # Build four translated blocks from current_T
    blocks = []
    for s in starts:
      base = s * delta - lo
      block = [(l + base, r + base) for (l, r) in current_T]
      blocks.append(block)

    # Assemble S with optional interleaving and reverse order to mix colors
    S = []
    if do_interleave:
      maxlen = max(len(b) for b in blocks)
      order = list(range(4))
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

    # Classic four connectors (caps)
    s0, s1, s2, s3 = starts
    S.append(((s0 - 1) * delta, (s1 - 1) * delta))  # left cap
    S.append(((s2 + 2) * delta, (s3 + 2) * delta))  # right cap
    S.append(((s0 + 2) * delta, (s2 - 1) * delta))  # cross 1
    S.append(((s1 + 2) * delta, (s3 - 1) * delta))  # cross 2
    return S

  def _thin_seed_with_offset(T, target_sz, offset):
    n = len(T)
    if n == 0 or target_sz <= 0:
      return []
    step = max(1, n // target_sz)
    off = offset % step if step > 0 else 0
    U = []
    i = off
    while len(U) < target_sz and i < n:
      U.append(T[i])
      i += step
    return U

  def _micro_round(current_T, windows, seed_tag=0, reverse_blocks=False, budget=0):
    """
    Take a thin seed U and translate it to each window; interleave blocks; add fractional connectors.
    Windows are expressed as fractions of the global span [glo, ghi].
    """
    if not current_T or budget <= 8:
      return []

    glo, ghi, G = _span_delta(current_T)

    # Seed size capped to keep omega modest; deterministic offset for reproducibility
    target_sz = max(8, min(40, len(current_T) // 260))
    BASE = 1337
    M = 2654435761  # Knuth multiplicative hash
    hashed = (BASE ^ (seed_tag * 374761393) ^ M) & 0xFFFFFFFF
    U = _thin_seed_with_offset(current_T, target_sz, offset=hashed)
    if not U:
      return []

    ulo = min(l for l, _ in U)

    # Build micro-blocks aligned to window left endpoints
    blocks = []
    for (fa, fb) in windows:
      a = max(0.05, min(0.90, fa))
      b = max(0.10, min(0.95, fb))
      if b <= a:
        b = a + 0.05
      win_lo = glo + int(round(a * G))
      base = win_lo - ulo
      block = [(l + base, r + base) for (l, r) in U]
      blocks.append(block)

    # Interleave micro-blocks with optional reversed block order
    order = list(range(len(blocks)))
    if reverse_blocks:
      order.reverse()
    interleaved = []
    maxlen = max(len(b) for b in blocks)
    for i in range(maxlen):
      for idx in order:
        b = blocks[idx]
        if i < len(b):
          interleaved.append(b[i])

    # Deterministic fractional-span connectors, including long-range cross4 ties
    con = [
      (glo + int(round(0.08 * G)), glo + int(round(0.30 * G))),  # left cap
      (glo + int(round(0.60 * G)), glo + int(round(0.92 * G))),  # right cap
      (glo + int(round(0.26 * G)), glo + int(round(0.56 * G))),  # cross 1
      (glo + int(round(0.44 * G)), glo + int(round(0.78 * G))),  # cross 2
      (glo + int(round(0.18 * G)), glo + int(round(0.84 * G))),  # cross4a
      (glo + int(round(0.12 * G)), glo + int(round(0.88 * G))),  # cross4b
    ]
    for a, b in con:
      if b > a:
        interleaved.append((a, b))

    if len(interleaved) > budget:
      interleaved = interleaved[:budget]
    return interleaved

  # Seed with one unit interval. We deliberately cap the spine depth to 5 (not 6) to
  # reduce omega while retaining strong FF pressure, then rely on two micro-phases.
  T = [(0, 1)]
  spine_depth = min(5, max(0, int(rounds)))  # guard user-provided rounds

  for ridx in range(spine_depth):
    # Alternate templates deterministically
    starts = template_bank[ridx % len(template_bank)]
    # Use sequential blocks; reverse order on odd rounds for mixing without big omega spikes
    do_interleave = False
    reverse_order = (ridx % 2 == 1)
    # Predict capacity
    if 4 * len(T) + 4 > CAP:
      break
    T = _apply_round(T, starts, do_interleave=do_interleave, reverse_order=reverse_order)

  # If close to capacity, skip micro-phases
  if len(T) >= CAP - 16:
    return T

  # Micro-phase A: interior windows; forward order
  windows_A = [(0.12, 0.22), (0.35, 0.45), (0.58, 0.68), (0.80, 0.90)]
  room = CAP - len(T)
  if room > 12:
    budgetA = min(room // 2, 300)  # conservative cap
    microA = _micro_round(T, windows_A, seed_tag=1, reverse_blocks=False, budget=budgetA)
    if microA:
      add = min(len(microA), CAP - len(T))
      T.extend(microA[:add])

  # Micro-phase B: distinct windows; reversed block order; includes long-range ties
  room = CAP - len(T)
  if room > 12:
    windows_B = [(0.05, 0.15), (0.28, 0.38), (0.60, 0.70), (0.82, 0.92)]
    budgetB = min(room, 350)  # fill further but still conservative
    microB = _micro_round(T, windows_B, seed_tag=2, reverse_blocks=True, budget=budgetB)
    if microB:
      add = min(len(microB), CAP - len(T))
      T.extend(microB[:add])

  # Optional extra micro iterations if requested; clamp to at most two more light passes
  extra_steps = max(0, int(phase2_iters) - 2)
  extra_steps = min(extra_steps, 2)
  for k in range(extra_steps):
    room = CAP - len(T)
    if room <= 12:
      break
    # Alternate window families each extra step
    if k % 2 == 0:
      micro = _micro_round(T, windows_A, seed_tag=3 + k, reverse_blocks=(k % 2 == 1), budget=min(room, 200))
    else:
      micro = _micro_round(T, windows_B, seed_tag=3 + k, reverse_blocks=(k % 2 == 1), budget=min(room, 200))
    if not micro:
      break
    T.extend(micro[: CAP - len(T)])

  # Hard cap
  if len(T) > CAP:
    T = T[:CAP]
  return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()