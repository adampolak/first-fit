# EVOLVE-BLOCK-START

def construct_intervals(rounds=6,
                        rotate_starts=True,
                        reverse_block_parity=True,
                        interleave_blocks=True,
                        phase2_iters=1):
  """
  Deterministic KT-style spine with a rotating-template scaffold and a thin micro-phase.

  Parameters:
    rounds (int): main expansion depth; near 6 yields ~9556 intervals for a single-seed KT spine.
    rotate_starts (bool): rotate among strong start templates when True.
    reverse_block_parity (bool): if True, flips even/odd interleaving parity each round.
    interleave_blocks (bool): enable interleaving on selected rounds.
    phase2_iters (int): micro-round iterations (safeguarded to at most two and budget-limited).

  Returns:
    intervals: list of (l, r) integer tuples, open intervals, in FF presentation order.
  """

  # Hard capacity guard to keep the total count < 10000
  CAP = 9800

  # Four strong start-pattern templates (rotated across rounds).
  spine_starts = (2, 6, 10, 14)
  template_bank = [
    (2, 6, 10, 14),  # classic KT
    (1, 5, 9, 13),   # left-shifted
    (3, 7, 11, 15),  # right-shifted
    (4, 8, 12, 16),  # stretched-right
  ]

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
    # Optional long-range cross connector to couple distant blocks (final rounds only)
    if add_cross4:
      S.append(((s0 + 4) * delta, (s3 + 4) * delta))

  def apply_round(current_T, starts, do_interleave=False, reverse_order=False, add_cross4=False):
    lo, hi, delta = _span_delta(current_T)

    # Build four translated blocks
    blocks = []
    for s in starts:
      base = s * delta - lo
      block = [(l + base, r + base) for (l, r) in current_T]
      blocks.append(block)

    # Build S either interleaving or sequential
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

    # Append classic connectors (+ optional long-range cross4)
    _append_connectors(S, starts, delta, add_cross4=add_cross4)
    return S

  # Stage 1: KT spine with rotating templates and parity-based interleaving.
  depth, _ = max_rounds_within_cap(len(T), rounds)
  for ridx in range(depth):
    starts = template_bank[ridx % len(template_bank)] if rotate_starts else spine_starts
    # Even/odd interleaving policy; optionally reverse block order each round.
    do_inter = bool(interleave_blocks and (ridx % 2 == 0))
    rev = bool(reverse_block_parity and (ridx % 2 == 1))
    # Enable a long-range cross connector only in the final two rounds to intensify color mixing
    add_c4 = (ridx >= depth - 2)
    T = apply_round(T, starts, do_interleave=do_inter, reverse_order=rev, add_cross4=add_c4)

  # If we are close to capacity, return the strong baseline.
  if len(T) >= CAP - 8:
    return T

  # Micro-phase A: insert a tiny tail of long caps near the end to boost FF mixing.
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
  caps = [cap_at(0.08, 0.60), cap_at(0.25, 0.75), cap_at(0.75, 0.92)]
  # Capacity-guarded insertion
  room = CAP - len(T)
  if room > 0:
    T = _insert_near_tail(T, caps[:room])

  if len(T) >= CAP - 16:
    return T

  # Stage 2: delta2-driven micro rounds using thin evenly-spaced seeds.
  def build_micro_delta_round(current_T, budget, iter_id=0, scheme=0):
    if not current_T or budget <= 8:
      return []

    glo = min(l for l, r in current_T)
    ghi = max(r for l, r in current_T)
    G = max(1, ghi - glo)

    # Thin seed from current_T, evenly spaced.
    seed_sz = max(8, min(32, len(current_T) // 300))
    stride = max(1, len(current_T) // max(1, seed_sz))
    U = [current_T[i] for i in range(0, len(current_T), stride)][:seed_sz]
    if not U:
      return []

    ulo = min(l for l, r in U)

    # Choose window profile
    if scheme == 1:
      # Distinct late-stage window set to diversify FF pressure
      window_fracs = [(0.05, 0.15), (0.28, 0.38), (0.60, 0.70), (0.82, 0.92)]
    else:
      # Slightly shift the window positions across iterations for diversification.
      shift = (iter_id % 3) * 0.02
      window_fracs = [
        (0.12 + shift, 0.22 + shift),
        (0.35 + shift, 0.45 + shift),
        (0.58 + shift, 0.68 + shift),
        (0.80 + shift, 0.90 + shift),
      ]
    # Clamp windows inside (0.05, 0.95)
    window_fracs = [
      (max(0.05, min(0.90, a)), max(0.10, min(0.95, b)))
      for (a, b) in window_fracs
    ]

    # Build translated micro-blocks aligned to these windows
    blocks = []
    for (fa, fb) in window_fracs:
      win_lo = glo + int(round(fa * G))
      base = win_lo - ulo
      block = [(l + base, r + base) for (l, r) in U]
      # Alternate internal reversal by block to break symmetry
      if ((int(round(fa * 100)) // 5) + iter_id + scheme) % 2 == 1:
        block = list(reversed(block))
      blocks.append(block)

    # Interleave micro-blocks (forward on even iter, reverse on odd iter)
    micro = []
    maxlen = max(len(b) for b in blocks)
    block_order = list(range(len(blocks)))
    if iter_id % 2 == 1:
      block_order.reverse()
    for i in range(maxlen):
      for idx in block_order:
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
    # Add a single long-range cross for the distinct scheme to enhance coupling
    if scheme == 1:
      micro_connectors.append((glo + int(round(0.12 * G)), glo + int(round(0.88 * G))))
    for a, b in micro_connectors:
      if b > a:
        micro.append((a, b))

    # Trim to available budget
    if len(micro) > budget:
      micro = micro[:budget]
    return micro

  # Execute up to two micro-rounds, capacity-guarded
  steps = min(max(0, int(phase2_iters)), 2)
  for iter_id in range(steps):
    room = CAP - len(T)
    if room <= 8:
      break
    micro = build_micro_delta_round(T, room, iter_id=iter_id)
    if not micro:
      break
    # Capacity guard and append
    avail = CAP - len(T)
    if len(micro) > avail:
      micro = micro[:avail]
    T.extend(micro)

  # Additional guarded micro-phase with a distinct window profile (scheme=1)
  room = CAP - len(T)
  if room > 8:
    micro2 = build_micro_delta_round(T, room, iter_id=steps, scheme=1)
    if micro2:
      avail = CAP - len(T)
      if len(micro2) > avail:
        micro2 = micro2[:avail]
      T.extend(micro2)

  return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()