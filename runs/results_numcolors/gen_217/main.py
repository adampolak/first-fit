# EVOLVE-BLOCK-START

def construct_intervals(rounds=6,
                        rotate_starts=True,
                        reverse_block_parity=True,
                        interleave_blocks=True,
                        phase2_iters=1):
  """
  Deterministic KT-style spine with dual micro-phases and a long-range cross4 layer.

  Parameters (kept for interface compatibility):
    rounds (int): main expansion depth; near 6 yields ~9556 intervals for a single-seed KT spine.
    rotate_starts (bool): rotate among strong start templates when True.
    reverse_block_parity (bool): flips even/odd interleaving parity each round when True.
    interleave_blocks (bool): enable interleaving on selected rounds when True.
    phase2_iters (int): number of primary micro-phase iterations (safeguarded to at most two).

  Returns:
    intervals: list of (l, r) integer tuples, open intervals, in FF presentation order.
  """

  # Hard capacity guard to keep the total count < 10000
  CAP = 9800

  # Global deterministic base seed for parity/window jitter derivation
  BASE_SEED = 91138241

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

  # Deterministic per-round seed derivation
  def derive_seed(base, stage_id, tag):
    x = (base ^ (stage_id * 0x9E3779B1) ^ (tag * 0x85EBCA6B)) & 0xFFFFFFFF
    # Extract small deterministic jitters/parities
    jitter = ((x >> 8) & 31) / 1000.0   # in [0, 0.031]
    parity = (x >> 3) & 1
    idx = (x % len(template_bank))
    return jitter, parity, idx

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

  def _append_connectors(S, starts, delta, add_cross4=False, round_parity=0):
    # Classic four connectors; preserves strong FF pressure while keeping omega modest.
    s0, s1, s2, s3 = starts
    S.append(((s0 - 1) * delta, (s1 - 1) * delta))  # left cap
    S.append(((s2 + 2) * delta, (s3 + 2) * delta))  # right cap
    S.append(((s0 + 2) * delta, (s2 - 1) * delta))  # cross 1
    S.append(((s1 + 2) * delta, (s3 - 1) * delta))  # cross 2
    # Optional long cross4 connector gated by parity to avoid systematic cliques
    if add_cross4 and (round_parity % 2 == 1):
      S.append(((s0 + 4) * delta, (s3 + 4) * delta))

  def apply_round(current_T, starts, do_interleave=False, reverse_order=False, ridx=0):
    lo, hi, delta = _span_delta(current_T)

    # Build four translated blocks
    blocks = []
    for s in starts:
      base = s * delta - lo
      block_src = current_T
      blocks.append([(l + base, r + base) for (l, r) in block_src])

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

    # Append classic connectors (plus a parity-gated long cross4 on odd rounds)
    _, parity_flag, _ = derive_seed(BASE_SEED, ridx, tag=1)
    _append_connectors(S, starts, delta, add_cross4=True, round_parity=parity_flag)
    return S

  # Stage 1: KT spine with rotating templates and parity-based interleaving.
  depth, _ = max_rounds_within_cap(len(T), rounds)
  for ridx in range(depth):
    _, parity_flag, idx = derive_seed(BASE_SEED, ridx, tag=0)
    starts = template_bank[idx] if rotate_starts else spine_starts
    # Even/odd interleaving policy; optionally reverse block order each round (both parity-gated)
    do_inter = bool(interleave_blocks and ((ridx + parity_flag) % 2 == 0))
    rev = bool(reverse_block_parity and ((ridx + parity_flag) % 2 == 1))
    T = apply_round(T, starts, do_interleave=do_inter, reverse_order=rev, ridx=ridx)

  # If we are close to capacity, return the strong baseline.
  if len(T) >= CAP - 24:
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

  if len(T) >= CAP - 64:
    return T

  # Stage 2: delta2-driven micro rounds using thin evenly-spaced seeds.
  def build_micro_delta_round(current_T, budget, iter_id=0, alt=False):
    if not current_T or budget <= 8:
      return []

    glo = min(l for l, r in current_T)
    ghi = max(r for l, r in current_T)
    G = max(1, ghi - glo)

    # Thin seed from current_T, evenly spaced (slightly larger for alt to increase coupling).
    base_seed = max(8, min(40, len(current_T) // 300))
    seed_sz = base_seed if not alt else max(10, min(48, len(current_T) // 250))
    stride = max(1, len(current_T) // max(1, seed_sz))
    U = [current_T[i] for i in range(0, len(current_T), stride)][:seed_sz]
    if not U:
      return []

    ulo = min(l for l, r in U)

    # Two families of windows: primary (A) with tiny deterministic jitter; alternate (B) fixed.
    jitter, parity_flag, _ = derive_seed(BASE_SEED, iter_id, tag=2 if not alt else 3)
    if not alt:
      shift = ((iter_id % 3) * 0.02) + 0.5 * jitter
      window_fracs = [
        (0.12 + shift, 0.22 + shift),
        (0.35 + shift, 0.45 + shift),
        (0.58 + shift, 0.68 + shift),
        (0.80 + shift, 0.90 + shift),
      ]
      window_fracs = [
        (max(0.05, min(0.90, a)), max(0.10, min(0.95, b)))
        for (a, b) in window_fracs
      ]
    else:
      window_fracs = [(0.05, 0.15), (0.28, 0.38), (0.60, 0.70), (0.82, 0.92)]

    # Build translated micro-blocks aligned to these windows
    blocks = []
    for (fa, fb) in window_fracs:
      win_lo = glo + int(round(fa * G))
      base = win_lo - ulo
      block = [(l + base, r + base) for (l, r) in U]
      # Alternate internal reversal by block to break symmetry
      if ((int(round(fa * 100)) // 5) + iter_id + parity_flag) % 2 == 1:
        block = list(reversed(block))
      blocks.append(block)

    # Interleave micro-blocks (forward on even tag, reverse on odd tag)
    micro = []
    maxlen = max(len(b) for b in blocks)
    block_order = list(range(len(blocks)))
    if (iter_id + parity_flag) % 2 == 1:
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
    # Add a longer-range cross4 only for the alternate micro-phase
    if alt:
      micro_connectors.append((glo + int(round(0.18 * G)), glo + int(round(0.84 * G))))
    for a, b in micro_connectors:
      if b > a:
        micro.append((a, b))

    # Trim to available budget
    if len(micro) > budget:
      micro = micro[:budget]
    return micro

  # Execute primary micro-phase (A) with up to two rounds, capacity-guarded
  steps = min(max(0, int(phase2_iters)), 2)
  for iter_id in range(steps):
    room = CAP - len(T)
    if room <= 12:
      break
    micro = build_micro_delta_round(T, room, iter_id=iter_id, alt=False)
    if not micro:
      break
    avail = CAP - len(T)
    if len(micro) > avail:
      micro = micro[:avail]
    T.extend(micro)

  # Secondary micro-phase (B) with distinct fixed windows, single guarded pass
  room = CAP - len(T)
  if room > 12:
    microB = build_micro_delta_round(T, room, iter_id=steps, alt=True)
    if microB:
      avail = CAP - len(T)
      if len(microB) > avail:
        microB = microB[:avail]
      T.extend(microB)

  # Long-range cross4 connector layer (thin), added last and strictly CAP-gated
  # These tie distant parts of the global span to force late FF color openings.
  if T and (CAP - len(T) >= 4):
    glo = min(l for l, r in T)
    ghi = max(r for l, r in T)
    G = max(1, ghi - glo)
    cross4_layer = [
      (glo + int(round(0.06 * G)), glo + int(round(0.88 * G))),
      (glo + int(round(0.14 * G)), glo + int(round(0.76 * G))),
      (glo + int(round(0.22 * G)), glo + int(round(0.64 * G))),
      (glo + int(round(0.30 * G)), glo + int(round(0.94 * G))),
    ]
    # Ensure positive length and do not exceed CAP
    extra = []
    for a, b in cross4_layer:
      if b > a:
        extra.append((a, b))
      if len(T) + len(extra) >= CAP:
        break
    T.extend(extra)

  # Final size cap
  if len(T) > CAP:
    T = T[:CAP]

  return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()