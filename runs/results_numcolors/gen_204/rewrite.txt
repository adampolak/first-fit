# EVOLVE-BLOCK-START

def construct_intervals(rounds=6,
                        rotate_starts=True,
                        reverse_block_parity=True,
                        interleave_blocks=True,
                        phase2_iters=1):
  """
  Deterministic KT-style spine with a rotating-template scaffold and a dual micro-phase.

  Parameters:
    rounds (int): main expansion depth; near 6 yields ~9556 intervals for a single-seed KT spine.
    rotate_starts (bool): rotate among strong start templates when True.
    reverse_block_parity (bool): if True, flips even/odd interleaving parity each round.
    interleave_blocks (bool): enable interleaving on selected rounds.
    phase2_iters (int): micro-round iterations (safeguarded and budget-limited; primary phase only).

  Returns:
    intervals: list of (l, r) integer tuples, open intervals, in FF presentation order.
  """

  # Hard capacity guard to keep the total count < 10000
  CAP = 9800

  # Four strong start-pattern templates (rotated across rounds).
  SPINE_STARTS = (2, 6, 10, 14)
  TEMPLATE_BANK = [
    (2, 6, 10, 14),  # classic KT
    (1, 5, 9, 13),   # left-shifted
    (3, 7, 11, 15),  # right-shifted
    (4, 8, 12, 16),  # stretched-right
  ]

  # Single global seed for deterministic tiny shifts in micro windows
  BASE_SEED = 91138217

  # Seed with one unit interval to allow six KT rounds within CAP.
  T = [(0, 1)]

  # Predictive size accounting to cap the number of full rounds.
  # KT growth per round: size -> 4*size + 4
  def round_next_size(sz):
    return 4 * sz + 4

  def max_rounds_within_cap(initial_size, max_rounds):
    sz = initial_size
    done = 0
    for _ in range(max(0, int(max_rounds))):
      nxt = round_next_size(sz)
      if nxt > CAP:
        break
      sz = nxt
      done += 1
    return done, sz

  def _span_delta(current_T):
    lo = min(l for l, r in current_T)
    hi = max(r for l, r in current_T)
    delta = hi - lo
    if delta <= 0:
      delta = 1
    return lo, hi, delta

  def _append_connectors(S, starts, delta, add_cross4=False):
    # Classic four connectors; preserves strong FF pressure while keeping omega modest.
    s0, s1, s2, s3 = starts
    S.append(((s0 - 1) * delta, (s1 - 1) * delta))  # left cap
    S.append(((s2 + 2) * delta, (s3 + 2) * delta))  # right cap
    S.append(((s0 + 2) * delta, (s2 - 1) * delta))  # cross 1
    S.append(((s1 + 2) * delta, (s3 - 1) * delta))  # cross 2
    if add_cross4:
      # A single long-range cross that ties distant blocks; gated to avoid raising omega much
      S.append(((s0 + 4) * delta, (s3 + 4) * delta))

  def apply_round(current_T, starts, do_interleave=False, reverse_order=False, add_cross4=False):
    lo, hi, delta = _span_delta(current_T)

    # Build four translated blocks; optionally reverse odd blocks to break symmetry
    blocks = []
    for b_idx, s in enumerate(starts):
      base = s * delta - lo
      block_src = current_T[::-1] if (reverse_order and (b_idx % 2 == 1)) else current_T
      block = [(l + base, r + base) for (l, r) in block_src]
      blocks.append(block)

    # Build S either interleaving or sequential
    S = []
    if do_interleave:
      maxlen = max(len(b) for b in blocks)
      order = list(range(4))
      # Parity-driven reversal of block order
      if reverse_order:
        order.reverse()
      # Apply a small rotation to diversify coupling
      rot = (len(current_T) + len(blocks[0])) % 4
      order = order[rot:] + order[:rot]
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

    # Append connectors; optional cross4 on alternate rounds
    _append_connectors(S, starts, delta, add_cross4=add_cross4)
    return S

  # Stage 1: KT spine with rotating templates and parity-based interleaving.
  depth, _ = max_rounds_within_cap(len(T), rounds)
  for ridx in range(depth):
    starts = TEMPLATE_BANK[ridx % len(TEMPLATE_BANK)] if rotate_starts else SPINE_STARTS
    # Even/odd interleaving policy; optionally reverse block order each round.
    do_inter = bool(interleave_blocks and (ridx % 2 == 0))
    rev = bool(reverse_block_parity and (ridx % 2 == 1))
    # Add a single cross4 connector layer on odd rounds (capacity-safe)
    add_cross4 = (ridx % 2 == 1)
    T = apply_round(T, starts, do_interleave=do_inter, reverse_order=rev, add_cross4=add_cross4)

  # If we are very close to capacity, return the strong baseline.
  if len(T) >= CAP - 8:
    return T

  # Micro-tail: insert a tiny set of long caps near the end to boost FF mixing.
  def _insert_near_tail(seq, intervals):
    out = list(seq)
    for i, iv in enumerate(intervals):
      pos = len(out) - (i * 2 + 1)
      if pos < 0:
        out.append(iv)
      else:
        out.insert(pos, iv)
    return out

  lo, hi, delta = _span_delta(T)
  # Long caps positioned as fractions of the current span; ensure monotone endpoints.
  def cap_at(a_frac, b_frac):
    L = lo + max(1, int(round(a_frac * delta)))
    R = lo + max(1, int(round(b_frac * delta)))
    if R <= L:
      R = L + 1
    return (L, R)
  tail_caps = [cap_at(0.08, 0.60), cap_at(0.25, 0.75), cap_at(0.75, 0.92)]
  room = CAP - len(T)
  if room > 0:
    T = _insert_near_tail(T, tail_caps[:min(len(tail_caps), room)])

  if len(T) >= CAP - 24:
    return T

  # Stage 2: dual micro phases using thin evenly-spaced seeds.
  def _thin_seed(current_T, max_seed):
    """Evenly spaced sample of current_T of size <= max_seed (deterministic)."""
    n = len(current_T)
    if n == 0 or max_seed <= 0:
      return []
    step = max(1, n // max_seed)
    return current_T[::step][:max_seed]

  def _micro_windows(iter_id, glo, G, alt=False):
    """Primary or alternate window families with tiny deterministic shifts."""
    if not alt:
      # Derive a tiny shift from BASE_SEED and iter_id; keep within [0, 0.02]
      shift_units = ((BASE_SEED ^ (iter_id * 1315423911)) % 3)
      shift = 0.01 * shift_units
      fracs = [
        (0.12 + shift, 0.22 + shift),
        (0.35 + shift, 0.45 + shift),
        (0.58 + shift, 0.68 + shift),
        (0.80 + shift, 0.90 + shift),
      ]
      # Clamp to keep safe margins
      fracs = [(max(0.05, a), min(0.95, b)) for (a, b) in fracs]
    else:
      # Distinct window set as recommended; fixed to stabilize behavior
      fracs = [(0.05, 0.15), (0.28, 0.38), (0.60, 0.70), (0.82, 0.92)]
    return fracs

  def _build_micro_delta_round(current_T, budget, iter_id=0, alt=False):
    if not current_T or budget <= 8:
      return []

    glo = min(l for l, r in current_T)
    ghi = max(r for l, r in current_T)
    G = max(1, ghi - glo)

    # Thin seed from current_T, evenly spaced (keep micro blocks small)
    seed_sz = max(8, min(48, len(current_T) // 250))
    U = _thin_seed(current_T, seed_sz)
    if not U:
      return []

    ulo = min(l for l, r in U)

    # Build translated micro-blocks aligned to windows
    window_fracs = _micro_windows(iter_id, glo, G, alt=alt)
    blocks = []
    for (fa, fb) in window_fracs:
      win_lo = glo + int(round(fa * G))
      base = win_lo - ulo
      block = [(l + base, r + base) for (l, r) in U]
      # Deterministic internal reversal to break symmetry
      tag = iter_id + (1 if alt else 0)
      if ((int(round(fa * 100)) // 5) + tag) % 2 == 1:
        block = list(reversed(block))
      blocks.append(block)

    # Interleave micro-blocks; reverse order on odd tag
    micro = []
    maxlen = max(len(b) for b in blocks)
    order = list(range(len(blocks)))
    tag = iter_id + (1 if alt else 0)
    if tag % 2 == 1:
      order.reverse()
    for i in range(maxlen):
      for idx in order:
        blk = blocks[idx]
        if i < len(blk):
          micro.append(blk[i])

    # Deterministic connectors across windows (fractional-span analog of KT caps)
    micro_connectors = [
      (glo + int(round(0.08 * G)), glo + int(round(0.30 * G))),  # left cap
      (glo + int(round(0.60 * G)), glo + int(round(0.92 * G))),  # right cap
      (glo + int(round(0.26 * G)), glo + int(round(0.56 * G))),  # cross 1
      (glo + int(round(0.44 * G)), glo + int(round(0.78 * G))),  # cross 2
    ]
    # A longer-range cross only in the alternate micro-phase
    if alt:
      micro_connectors.append((glo + int(round(0.18 * G)), glo + int(round(0.84 * G))))

    for a, b in micro_connectors:
      if b > a:
        micro.append((a, b))

    # Trim to available budget
    if len(micro) > budget:
      micro = micro[:budget]
    return micro

  # Primary micro-phase: up to two iterations (guarded by phase2_iters and CAP)
  steps = min(max(0, int(phase2_iters)), 2)
  for iter_id in range(steps):
    room = CAP - len(T)
    if room <= 8:
      break
    micro = _build_micro_delta_round(T, room, iter_id=iter_id, alt=False)
    if not micro:
      break
    avail = CAP - len(T)
    if len(micro) > avail:
      micro = micro[:avail]
    T.extend(micro)

  # Secondary micro-phase with distinct windows to capture missed interactions (one pass)
  room = CAP - len(T)
  if room > 8:
    microB = _build_micro_delta_round(T, room, iter_id=steps, alt=True)
    if microB:
      avail = CAP - len(T)
      if len(microB) > avail:
        microB = microB[:avail]
      T.extend(microB)

  # Long-range global cross4 connector layer (very thin; add at most 3)
  # Ties distant positions at fractional offsets across the entire span.
  room = CAP - len(T)
  if room > 0 and T:
    glo = min(l for l, r in T)
    ghi = max(r for l, r in T)
    G = max(1, ghi - glo)
    cross4s = [
      (glo + int(round(0.07 * G)), glo + int(round(0.91 * G))),
      (glo + int(round(0.19 * G)), glo + int(round(0.83 * G))),
      (glo + int(round(0.31 * G)), glo + int(round(0.77 * G))),
    ]
    # Append up to 'room' connectors
    for iv in cross4s[:room]:
      a, b = iv
      if b > a:
        T.append(iv)
        room -= 1
        if room == 0:
          break

  return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()