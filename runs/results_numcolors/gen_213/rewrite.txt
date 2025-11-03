# EVOLVE-BLOCK-START

def construct_intervals(seed_count=1):
  """
  Construct a sequence of open intervals (l, r) in FirstFit arrival order.
  Same interface as the original: only seed_count is accepted (and safely handled).

  Returns:
    intervals: list of (l, r) integer tuples.
  """
  CAP = 9800

  # Deterministic KT spine templates (rotated to diversify coupling)
  template_bank = [
    (2, 6, 10, 14),  # classic KT
    (1, 5, 9, 13),   # left-shifted
    (3, 7, 11, 15),  # right-shifted
    (4, 8, 12, 16),  # stretched-right
  ]

  # Seed: single interval to keep early omega minimal
  if seed_count <= 1:
    T = [(0, 1)]
  else:
    # Deterministic sparse multi-seed; capped to avoid early omega blow-up
    step = 3
    seeds = max(1, int(seed_count))
    T = [(i * step, i * step + 1) for i in range(seeds)]
    T = T[:4]

  def _span_delta(current_T):
    lo = min(l for l, r in current_T)
    hi = max(r for l, r in current_T)
    delta = hi - lo
    if delta <= 0:
      delta = 1
    return lo, hi, delta

  def _append_connectors(S, starts, delta):
    # Classic KT connectors (keep omega modest)
    s0, s1, s2, s3 = starts
    S.append(((s0 - 1) * delta, (s1 - 1) * delta))  # left cap
    S.append(((s2 + 2) * delta, (s3 + 2) * delta))  # right cap
    S.append(((s0 + 2) * delta, (s2 - 1) * delta))  # cross 1
    S.append(((s1 + 2) * delta, (s3 - 1) * delta))  # cross 2

  def _apply_round(current_T, starts, do_interleave=False, reverse_order=False):
    lo, hi, delta = _span_delta(current_T)

    # Build four translated blocks
    blocks = []
    for s in starts:
      base = s * delta - lo
      block = [(l + base, r + base) for (l, r) in current_T]
      blocks.append(block)

    # Build S
    S = []
    if do_interleave:
      # Interleave with optional reversed block order
      order = list(range(4))
      if reverse_order:
        order.reverse()
      maxlen = max(len(b) for b in blocks)
      for i in range(maxlen):
        for idx in order:
          blk = blocks[idx]
          if i < len(blk):
            S.append(blk[i])
    else:
      # Sequential composition with optional reversed block order
      if reverse_order:
        blocks = list(reversed(blocks))
      for blk in blocks:
        S.extend(blk)

    _append_connectors(S, starts, delta)
    return S

  # Stage 1: Six KT spine rounds with rotating templates and parity-based interleaving
  rounds = 6
  for ridx in range(rounds):
    starts = template_bank[ridx % len(template_bank)]
    # Even rounds: interleave forward. Odd rounds: sequential with reversed block order.
    do_interleave = (ridx % 2 == 0)
    reverse_order = (ridx % 2 == 1)
    T = _apply_round(T, starts, do_interleave=do_interleave, reverse_order=reverse_order)
    if len(T) >= CAP:
      return T[:CAP]

  # Early return if at capacity
  if len(T) >= CAP - 8:
    return T[:CAP]

  # Stage 2: Deterministic micro-phase rounds
  # Two primary micro-rounds + one alternate micro-round (all CAP-guarded)

  def _thin_seed(current_T, max_seed):
    n = len(current_T)
    if n == 0 or max_seed <= 0:
      return []
    step = max(1, n // max_seed)
    U = current_T[::step][:max_seed]
    return U

  def _build_micro_delta_round(current_T, budget, iter_id=0, alt=False):
    if not current_T or budget <= 8:
      return []

    glo = min(l for l, r in current_T)
    ghi = max(r for l, r in current_T)
    G = max(1, ghi - glo)

    # Seed size tuned to be thin to avoid inflating omega
    base_seed = max(8, min(32, len(current_T) // 300))
    seed_sz = base_seed if not alt else max(8, min(40, len(current_T) // 250))
    U = _thin_seed(current_T, seed_sz)
    if not U:
      return []

    ulo = min(l for l, r in U)

    # Windows: primary A with slight shifts; alternate B fixed distinct windows
    if not alt:
      shift = (iter_id % 3) * 0.02
      window_fracs = [
        (0.12 + shift, 0.22 + shift),
        (0.35 + shift, 0.45 + shift),
        (0.58 + shift, 0.68 + shift),
        (0.80 + shift, 0.90 + shift),
      ]
      # Clamp inside (0.05, 0.95)
      window_fracs = [
        (max(0.05, min(0.90, a)), max(0.10, min(0.95, b)))
        for (a, b) in window_fracs
      ]
    else:
      window_fracs = [(0.05, 0.15), (0.28, 0.38), (0.60, 0.70), (0.82, 0.92)]

    # Build translated micro-blocks aligned to windows with parity-breaking reversals
    blocks = []
    for (fa, fb) in window_fracs:
      win_lo = glo + int(round(fa * G))
      base = win_lo - ulo
      block = [(l + base, r + base) for (l, r) in U]
      tag = iter_id if not alt else (iter_id + 1)
      if ((int(round(fa * 100)) // 5) + tag) % 2 == 1:
        block = list(reversed(block))
      blocks.append(block)

    # Interleave policy: forward on even tag, reverse on odd tag
    micro = []
    tag = iter_id if not alt else (iter_id + 1)
    block_order = list(range(len(blocks)))
    if tag % 2 == 1:
      block_order.reverse()
    maxlen = max(len(b) for b in blocks)
    for i in range(maxlen):
      for idx in block_order:
        blk = blocks[idx]
        if i < len(blk):
          micro.append(blk[i])

    # Deterministic connectors across windows (fractional-span analog of KT)
    micro_connectors = [
      (glo + int(round(0.08 * G)), glo + int(round(0.30 * G))),  # left cap
      (glo + int(round(0.60 * G)), glo + int(round(0.92 * G))),  # right cap
      (glo + int(round(0.26 * G)), glo + int(round(0.56 * G))),  # cross 1
      (glo + int(round(0.44 * G)), glo + int(round(0.78 * G))),  # cross 2
    ]
    if alt:
      # Long-range cross to strengthen long-distance coupling
      micro_connectors.append((glo + int(round(0.18 * G)), glo + int(round(0.84 * G))))
    for a, b in micro_connectors:
      if b > a:
        micro.append((a, b))

    # Trim to available budget
    if len(micro) > budget:
      micro = micro[:budget]
    return micro

  # Execute up to two primary micro-rounds
  for iter_id in range(2):
    room = CAP - len(T)
    if room <= 8:
      break
    micro = _build_micro_delta_round(T, room, iter_id=iter_id, alt=False)
    if not micro:
      break
    if len(micro) > room:
      micro = micro[:room]
    T.extend(micro)

  # Alternate micro-phase (distinct window family + long-range connector), CAP-guarded
  room = CAP - len(T)
  if room > 8:
    microB = _build_micro_delta_round(T, room, iter_id=2, alt=True)
    if microB:
      avail = CAP - len(T)
      if len(microB) > avail:
        microB = microB[:avail]
      T.extend(microB)

  # Final capacity guard
  if len(T) > CAP:
    T = T[:CAP]
  return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()