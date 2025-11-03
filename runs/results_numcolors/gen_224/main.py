# EVOLVE-BLOCK-START

def construct_intervals(seed_count=1):
  """
  Construct a sequence of open intervals for FirstFit in presentation order.
  Maintains the same I/O as prior versions: returns a list of (l, r) integer tuples.
  """

  # Hard capacity to keep intervals < 10000 and allow room for micro-phases
  CAP = 9800

  # Deterministic global seed to derive minor parity choices and offsets (no randomness)
  BASE_SEED = 137

  # Rotating start templates for the KT-style backbone (classic set)
  TEMPLATE_BANK = [
    (2, 6, 10, 14),  # T1: classic KT
    (1, 5, 9, 13),   # T2: left-shifted
    (3, 7, 11, 15),  # T3: right-shifted
    (4, 8, 12, 16),  # T4: stretched-right
  ]

  # Seed with a single unit interval to avoid early omega inflation
  T = [(0, 1)] if seed_count <= 1 else [(0, 1)]

  # Helpers
  def _span(TS):
    lo = min(l for l, r in TS)
    hi = max(r for l, r in TS)
    delta = hi - lo
    if delta <= 0:
      delta = 1
    return lo, hi, delta

  def _apply_round(current_T, starts, do_interleave=False, reverse_order=False, add_cross4=False):
    lo, hi, delta = _span(current_T)

    # Build four translated blocks
    blocks = []
    for s in starts:
      base = s * delta - lo
      block = [(l + base, r + base) for (l, r) in current_T]
      # Gentle internal reversal by start parity to break symmetry
      if ((s + BASE_SEED) // 2) % 2 == 1:
        block = list(reversed(block))
      blocks.append(block)

    # Assemble S: either interleave round-robin or keep sequential (optionally reversed)
    S = []
    if do_interleave:
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
      if reverse_order:
        blocks = list(reversed(blocks))
      for blk in blocks:
        S.extend(blk)

    # Classic four connectors
    s0, s1, s2, s3 = starts
    S.append(((s0 - 1) * delta, (s1 - 1) * delta))  # left cap
    S.append(((s2 + 2) * delta, (s3 + 2) * delta))  # right cap
    S.append(((s0 + 2) * delta, (s2 - 1) * delta))  # cross 1
    S.append(((s1 + 2) * delta, (s3 - 1) * delta))  # cross 2

    # Optional long-range connector to enhance cross-scale coupling
    if add_cross4:
      # Carefully offset to avoid coinciding with classic caps
      S.append(((s0 + 4) * delta, (s3 + 4) * delta))

    return S

  def _predict_next_size(sz):
    # Each round grows as 4*sz + 4 (four blocks + four connectors)
    return 4 * sz + 4

  def _thin_seed(seq, target):
    n = len(seq)
    if n == 0 or target <= 0:
      return []
    step = max(1, n // target)
    return seq[::step][:target]

  # Choose backbone depth to leave space for twin layer + micro phases
  # Leave at least ~512 slots for gadgets after backbone
  size = len(T)
  rounds = 0
  while rounds < 6:
    nxt = _predict_next_size(size)
    if nxt > CAP - 512:
      break
    size = nxt
    rounds += 1

  # Stage 1: Twin-spine backbone with parity-based interleaving and rotating templates
  for ridx in range(rounds):
    starts = TEMPLATE_BANK[ridx % len(TEMPLATE_BANK)]
    do_inter = (ridx % 2 == 0)
    rev = (ridx % 2 == 1)
    # Only add the long cross on the last executed round to limit clique growth
    add_c4 = (ridx == rounds - 1) and (rounds >= 3)
    T = _apply_round(T, starts, do_interleave=do_inter, reverse_order=rev, add_cross4=add_c4)
    if len(T) >= CAP:
      return T[:CAP]

    # Mid-round twin layer after rounds 2 and 4 (0-indexed), strictly CAP-gated
    if ridx in (1, 3) and len(T) + 128 < CAP:
      # This "twin" embeds a thin seed at half-scale offsets to couple colors
      glo, ghi, G = _span(T)
      delta2 = max(1, G // 2)
      seed_sz = max(8, min(48, len(T) // 320))
      U = _thin_seed(T, seed_sz)
      if U:
        ulo = min(l for l, r in U)
        starts_twin = TEMPLATE_BANK[(ridx + 1) % len(TEMPLATE_BANK)]
        blocks = []
        for s in starts_twin:
          base = s * delta2 - ulo + (BASE_SEED % 3)
          block = [(l + base, r + base) for (l, r) in U]
          if ((s + ridx) % 3) == 1:
            block = list(reversed(block))
          blocks.append(block)
        # Parity-based interleaving
        twin = []
        maxlen = max(len(b) for b in blocks)
        order = [0, 1, 2, 3]
        if (BASE_SEED + ridx) % 2 == 1:
          order.reverse()
        for i in range(maxlen):
          for idx in order:
            blk = blocks[idx]
            if i < len(blk):
              twin.append(blk[i])

        # Add sparse connectors at half scale to avoid omega blow-up
        s0, s1, s2, s3 = starts_twin
        connectors = [
          (glo + (s0 - 1) * delta2, glo + (s1 - 1) * delta2),
          (glo + (s2 + 2) * delta2, glo + (s3 + 2) * delta2),
        ]
        twin.extend(connectors)
        # Capacity guard
        room = CAP - len(T)
        if room > 0:
          twin = twin[:room]
          T.extend(twin)

  if len(T) >= CAP - 64:
    return T[:CAP]

  # Cross4 connector layer (deterministic, long-range, capped)
  # Ties distant positions ~4–6 delta apart to strengthen mixing across scales.
  lo, hi, delta = _span(T)
  span = max(1, hi - lo)
  cross4 = [
    (lo + 1 * delta, lo + 5 * delta),
    (lo + 2 * delta, lo + 6 * delta),
    (lo + 3 * delta, lo + 7 * delta),
  ]
  for c in cross4:
    if len(T) >= CAP:
      break
    a, b = c
    if b > a:
      T.append((a, b))

  if len(T) >= CAP - 64:
    return T[:CAP]

  # Micro-phase builder with two window families
  def _build_micro(current_T, budget, windows, iter_id=0):
    if not current_T or budget <= 8:
      return []

    glo = min(l for l, r in current_T)
    ghi = max(r for l, r in current_T)
    G = max(1, ghi - glo)
    # Thin, bounded seed
    seed_sz = max(8, min(40, len(current_T) // 250))
    U = _thin_seed(current_T, seed_sz)
    if not U:
      return []
    ulo = min(l for l, r in U)

    # Build translated micro-blocks aligned to fractional windows
    blocks = []
    for (fa, fb) in windows:
      win_lo = glo + int(round(fa * G))
      base = win_lo - ulo
      block = [(l + base, r + base) for (l, r) in U]
      # Internal reversal keyed by iter_id and window
      tag = ((int(fa * 100) + 3 * iter_id + BASE_SEED) // 5) % 2
      if tag == 1:
        block = list(reversed(block))
      blocks.append(block)

    # Interleave with parity on iter_id
    micro = []
    maxlen = max(len(b) for b in blocks)
    order = list(range(len(blocks)))
    if (iter_id + BASE_SEED) % 2 == 1:
      order.reverse()
    for i in range(maxlen):
      for idx in order:
        blk = blocks[idx]
        if i < len(blk):
          micro.append(blk[i])

    # Fractional-span connectors for the micro scale
    micro_connectors = [
      (glo + int(round(0.08 * G)), glo + int(round(0.30 * G))),  # left cap
      (glo + int(round(0.26 * G)), glo + int(round(0.56 * G))),  # cross 1
      (glo + int(round(0.44 * G)), glo + int(round(0.78 * G))),  # cross 2
      (glo + int(round(0.60 * G)), glo + int(round(0.92 * G))),  # right cap
    ]
    # Add a long cross only on odd iterations
    if iter_id % 2 == 1:
      micro_connectors.append((glo + int(round(0.18 * G)), glo + int(round(0.84 * G))))
    for a, b in micro_connectors:
      if b > a:
        micro.append((a, b))

    if len(micro) > budget:
      micro = micro[:budget]
    return micro

  # Two guarded micro-phases with distinct window sets
  W1 = [(0.12, 0.22), (0.35, 0.45), (0.58, 0.68), (0.80, 0.90)]
  W2 = [(0.05, 0.15), (0.28, 0.38), (0.60, 0.70), (0.82, 0.92)]

  # Micro-phase 1
  room = CAP - len(T)
  if room > 16:
    micro1 = _build_micro(T, room, W1, iter_id=0)
    if micro1:
      micro1 = micro1[:room]
      T.extend(micro1)

  # Gate Micro-phase 2 to stay within CAP and avoid omega spikes near capacity
  room = CAP - len(T)
  if room > 16:
    micro2 = _build_micro(T, room, W2, iter_id=1)
    if micro2:
      micro2 = micro2[:room]
      T.extend(micro2)

  # Final capacity trim and normalization
  if len(T) > CAP:
    T = T[:CAP]

  # Normalize to integer open intervals and ensure r > l
  out = []
  for (l, r) in T:
    li = int(l)
    ri = int(r)
    if ri <= li:
      ri = li + 1
    out.append((li, ri))

  return out

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()