# EVOLVE-BLOCK-START

def construct_intervals(rounds=6,
                        rotate_starts=True,
                        reverse_block_parity=True,
                        interleave_blocks=True,
                        phase2_iters=1):
  """
  Braided six-template KT-style spine with deterministic selection,
  post-spine connectors and densifier, and an adaptive two-pass micro-phase.

  Returns:
    list[(l, r)] integer open intervals in presentation order for FirstFit.
  """

  # Hard capacity guard
  CAP = 9800
  BASE_SEED = 2654435761  # fixed; deterministic

  # Six start-pattern templates; maintain 4-block expansions for omega control
  template_bank = [
    (2, 6, 10, 14),  # classic KT
    (1, 5, 9, 13),   # left-shifted
    (3, 7, 11, 15),  # right-shifted
    (4, 8, 12, 16),  # stretched-right
    (2, 5, 11, 15),  # asymmetric 1
    (1, 7, 10, 14),  # asymmetric 2
  ]
  default_starts = (2, 6, 10, 14)

  # Simple LCG for deterministic hashing
  def lcg_step(x):
    return (1103515245 * x + 12345) & 0xFFFFFFFF

  def round_hash(ridx):
    x = (BASE_SEED ^ (ridx + 1)) & 0xFFFFFFFF
    x = lcg_step(x)
    x = lcg_step(x)
    return x

  # Seed: single unit interval
  T = [(0, 1)]

  # Span helper
  def span_info(seq):
    lo = min(l for l, r in seq)
    hi = max(r for l, r in seq)
    delta = hi - lo
    if delta <= 0:
      delta = 1
    return lo, hi, delta

  # KT growth per round: n -> 4n + 4
  def next_size(sz):
    return 4 * sz + 4

  # Cap how many full rounds fit within CAP ignoring later phases
  def rounds_within_cap(init_sz, max_rounds):
    sz = init_sz
    done = 0
    for _ in range(max(0, int(max_rounds))):
      nxt = next_size(sz)
      # Allow headroom for post-spine additions
      if nxt > CAP - 128:
        break
      sz = nxt
      done += 1
    return done

  # Append classic four KT connectors
  def append_connectors(S, starts, delta):
    s0, s1, s2, s3 = starts
    S.append(((s0 - 1) * delta, (s1 - 1) * delta))  # left cap
    S.append(((s2 + 2) * delta, (s3 + 2) * delta))  # right cap
    S.append(((s0 + 2) * delta, (s2 - 1) * delta))  # cross 1
    S.append(((s1 + 2) * delta, (s3 - 1) * delta))  # cross 2

  # Apply one KT-style round with deterministic braid decisions
  def apply_round(current_T, ridx, starts, allow_interleave, allow_reverse):
    lo, hi, delta = span_info(current_T)

    # Build four translated copies
    blocks = []
    for s in starts:
      base = s * delta - lo
      block = [(l + base, r + base) for (l, r) in current_T]
      blocks.append(block)

    # Braid and ordering: hash plus parity control
    H = round_hash(ridx)
    do_inter = bool(allow_interleave and ((ridx % 2) == 0 or (H & 1) == 1))
    rev_order = bool(allow_reverse and ((ridx % 2) == 1 or (H & 2) == 2))

    # Reorder blocks deterministically
    order = list(range(4))
    if (H >> 3) & 1:
      order = [0, 2, 1, 3]
    if rev_order:
      order = list(reversed(order))

    S = []
    if do_inter:
      maxlen = max(len(b) for b in blocks)
      for i in range(maxlen):
        for idx in order:
          blk = blocks[idx]
          if i < len(blk):
            S.append(blk[i])
    else:
      for idx in order:
        S.extend(blocks[idx])

    # Append classic connectors
    append_connectors(S, starts, delta)
    return S

  # Stage 1: braided six-template spine
  depth = rounds_within_cap(len(T), rounds)
  for ridx in range(depth):
    if rotate_starts:
      H = round_hash(ridx)
      starts = template_bank[H % len(template_bank)]
    else:
      starts = default_starts
    T = apply_round(
      T, ridx, starts,
      allow_interleave=interleave_blocks,
      allow_reverse=reverse_block_parity
    )

  # Early exit if near capacity
  if len(T) >= CAP - 32:
    return T

  # Stage 1.5: deterministic long-range post-spine connectors (not replicated)
  lo, hi, delta = span_info(T)
  G = delta
  # Fractions chosen to avoid clustering; ensure integer endpoints and l<r
  lr_fracs = [
    (0.07, 0.31),
    (0.19, 0.83),
    (0.28, 0.52),
    (0.41, 0.73),
    (0.58, 0.94),
  ]
  post_caps = []
  for a, b in lr_fracs:
    L = lo + max(1, int(round(a * G)))
    R = lo + max(2, int(round(b * G)))
    if R > L:
      post_caps.append((L, R))
  room = CAP - len(T)
  if room > 0 and post_caps:
    if len(post_caps) > room:
      post_caps = post_caps[:room]
    T.extend(post_caps)

  if len(T) >= CAP - 32:
    return T

  # Stage 1.75: sparse densifier pins (very short, to avoid omega spikes)
  lo, hi, delta = span_info(T)
  G = delta
  dens_room = max(0, CAP - len(T))
  if dens_room > 0:
    # Evenly spaced sample from T
    sample_target = min(96, max(24, len(T) // 400))
    stride = max(1, len(T) // max(1, sample_target))
    U = [T[i] for i in range(0, len(T), stride)][:sample_target]
    # Short pins near midpoints, shifted into five window anchors
    eps = max(1, G // 1024)
    anchors = [0.12, 0.30, 0.50, 0.70, 0.88]
    pins = []
    for idx, (l, r) in enumerate(U):
      mid = (l + r) // 2
      frac = anchors[idx % len(anchors)]
      base = lo + int(round(frac * G))
      L = base + ((mid - base) // 8)  # pull toward base; keep short
      R = L + eps
      if R > L:
        pins.append((L, R))
    if pins:
      if len(pins) > dens_room:
        pins = pins[:dens_room]
      T.extend(pins)

  if len(T) >= CAP - 16:
    return T

  # Stage 2: adaptive, two-pass micro-phase
  glo, ghi, G = span_info(T)
  remaining = max(0, CAP - len(T))
  if remaining == 0:
    return T

  # Split budget across two passes (bounded by phase2_iters)
  pass_cnt = 2 if phase2_iters >= 1 else 1
  pass_cnt = min(pass_cnt, 2)
  budget_A = (remaining * 3) // 5  # ~60% to pass A
  budget_B = remaining - budget_A

  def micro_pass(current_T, budget, pass_id):
    if budget <= 8 or not current_T:
      return []
    glo, ghi, G = span_info(current_T)

    # Thin evenly-spaced seed
    base_seed = max(12, min(40, len(current_T) // 280))
    stride = max(1, len(current_T) // max(1, base_seed))
    U = [current_T[i] for i in range(0, len(current_T), stride)][:base_seed]
    if not U:
      return []

    ulo = min(l for l, r in U)

    # Five-window families; different families per pass
    if pass_id == 0:
      windows = [(0.10, 0.18), (0.28, 0.36), (0.46, 0.54), (0.64, 0.72), (0.82, 0.90)]
    else:
      windows = [(0.05, 0.12), (0.22, 0.30), (0.40, 0.48), (0.58, 0.66), (0.76, 0.88)]

    # Build translated micro-blocks; alternate internal order deterministically
    blocks = []
    for widx, (fa, fb) in enumerate(windows):
      win_lo = glo + int(round(fa * G))
      base = win_lo - ulo
      block = [(l + base, r + base) for (l, r) in U]
      if ((pass_id + widx) % 2) == 1:
        block.reverse()
      blocks.append(block)

    # Interleave with deterministic order to enhance mixing
    interleave_order = list(range(len(blocks)))
    if pass_id == 1:
      interleave_order = [0, 2, 1, 3, 4]
    micro = []
    maxlen = max(len(b) for b in blocks)
    for i in range(maxlen):
      for idx in interleave_order:
        blk = blocks[idx]
        if i < len(blk):
          micro.append(blk[i])

    # Fractional-span connectors within the micro-phase
    connectors = [
      (glo + int(round(0.08 * G)), glo + int(round(0.30 * G))),
      (glo + int(round(0.22 * G)), glo + int(round(0.56 * G))),
      (glo + int(round(0.44 * G)), glo + int(round(0.78 * G))),
      (glo + int(round(0.60 * G)), glo + int(round(0.92 * G))),
    ]
    for a, b in connectors:
      if b > a:
        micro.append((a, b))

    # Short edge pins (very limited) to avoid omega spikes
    eps = max(1, G // 2048)
    for frac in ([0.14, 0.86] if pass_id == 1 else [0.18, 0.82]):
      a = glo + int(round(frac * G))
      b = a + eps
      if b > a:
        micro.append((a, b))

    if len(micro) > budget:
      micro = micro[:budget]
    return micro

  # Execute passes capacity-guarded
  if pass_cnt >= 1 and budget_A > 0:
    microA = micro_pass(T, budget_A, pass_id=0)
    if microA:
      avail = CAP - len(T)
      if len(microA) > avail:
        microA = microA[:avail]
      T.extend(microA)

  if pass_cnt >= 2:
    avail = CAP - len(T)
    if avail > 0 and budget_B > 0:
      microB = micro_pass(T, min(avail, budget_B), pass_id=1)
      if microB:
        avail = CAP - len(T)
        if len(microB) > avail:
          microB = microB[:avail]
        T.extend(microB)

  return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()