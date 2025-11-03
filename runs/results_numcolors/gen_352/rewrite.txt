# EVOLVE-BLOCK-START

def construct_intervals(enable_alt_microphase=True):
  """
  Construct a deterministic sequence of open intervals for FirstFit.

  Same inputs/outputs as the prior program:
    enable_alt_microphase (bool): whether to run the second micro-phase pass.

  Returns:
    list of (l, r) integer tuples representing open intervals in arrival order.
  """

  # Capacity guard (keep < 10000)
  CAP = 9800

  # Expanded template bank (six templates) to rotate across rounds
  template_bank = [
    (2, 6, 10, 14),
    (1, 5, 9, 13),
    (3, 7, 11, 15),
    (4, 8, 12, 16),
    (5, 9, 13, 17),
    (6, 10, 14, 18),
  ]

  # Deterministic interleaving schedule: interleave on rounds with idx % 3 != 0,
  # reverse-block-order on rounds where idx % 4 == 1 (deterministic pattern).
  def should_interleave(round_idx):
    return (round_idx % 3) != 0

  def should_reverse(round_idx):
    return (round_idx % 4) == 1

  # Seed with a single unit interval
  T = [(0, 1)]

  # Utility: compute span and delta (nonzero)
  def span_info(intervals):
    lo = min(l for l, r in intervals)
    hi = max(r for l, r in intervals)
    delta = hi - lo
    if delta <= 0:
      delta = 1
    return lo, hi, delta

  # Apply a KT-like round with given starts (4-tuple)
  def apply_spine_round(current_T, starts, do_interleave=False, reverse_order=False):
    lo, hi, delta = span_info(current_T)

    # Build four translated blocks
    blocks = []
    for s in starts:
      base = s * delta - lo
      block = [(l + base, r + base) for (l, r) in current_T]
      blocks.append(block)

    # Determine block output order using deterministic ordering function
    S = []
    if do_interleave:
      order = sorted(range(4), key=lambda i: starts[i], reverse=reverse_order)
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

    # Classic four connectors (Figure-4 style) using delta-scaled positions
    s0, s1, s2, s3 = starts
    # left cap, right cap, cross1, cross2
    S.append(((s0 - 1) * delta, (s1 - 1) * delta))
    S.append(((s2 + 2) * delta, (s3 + 2) * delta))
    S.append(((s0 + 2) * delta, (s2 - 1) * delta))
    S.append(((s1 + 2) * delta, (s3 - 1) * delta))

    return S

  # Predictive round count: we don't want to overshoot CAP
  # Growth: size -> 4*size + 4 per round
  def rounds_we_can_do(initial_size, max_rounds):
    sz = initial_size
    done = 0
    for r in range(max_rounds):
      nxt = 4 * sz + 4
      if nxt > CAP:
        break
      sz = nxt
      done += 1
    return done

  # Stage 1: KT-like spine with rotation through 6 templates.
  max_spine_rounds = 8  # try up to 8, but will be clipped by CAP-aware check
  can_do = rounds_we_can_do(len(T), max_spine_rounds)
  for ridx in range(can_do):
    starts = template_bank[ridx % len(template_bank)]
    inter = should_interleave(ridx)
    rev = should_reverse(ridx)
    T = apply_spine_round(T, starts, do_interleave=inter, reverse_order=rev)

    # CAP-guard break
    if len(T) >= CAP - 12:
      break

    # Per-round deterministic density multiplier applied to a small central subset
    # (increases local pressure without exploding omega). Apply only on selected rounds.
    if ridx % 2 == 0 and len(T) < CAP - 32:
      # pick a small central window and clone-shrink a few intervals
      lo, hi, delta = span_info(T)
      center = lo + delta // 2
      # sample up to 6 intervals near the center
      sampled = [iv for iv in T if iv[0] <= center <= iv[1]]
      # fall back to uniformly spaced sample if none
      if not sampled:
        step = max(1, len(T) // 100)
        sampled = [T[i] for i in range(0, len(T), step)][:6]
      densified = []
      shrink = max(1, delta // 2048)
      for (l, r) in sampled[:6]:
        if r - l > 2 * shrink:
          densified.append((l + shrink, r - shrink))
      # insert densified intervals near the end (preserve arrival mixing)
      for iv in densified:
        if len(T) >= CAP:
          break
        T.append(iv)

  # Early exit if near capacity
  if len(T) >= CAP - 8:
    return T[:CAP]

  # Micro-phase budget splitting: deterministic two-pass split (approx 60/40)
  remaining = CAP - len(T)
  passA_budget = max(8, int(round(0.60 * remaining)))
  passB_budget = remaining - passA_budget

  # Micro-phase builder: thin-window replication plus local connectors.
  # Two families: primary (shifted) and alternate (offset windows + long cross).
  def build_micro_pass(current_T, budget, pass_id=0, alt=False):
    if not current_T or budget <= 8:
      return []

    glo = min(l for l, r in current_T)
    ghi = max(r for l, r in current_T)
    G = max(1, ghi - glo)

    # thin evenly spaced seed U
    base_seed = max(8, min(40, len(current_T) // (220 if alt else 300)))
    stride = max(1, len(current_T) // max(1, base_seed))
    U = [current_T[i] for i in range(0, len(current_T), stride)][:base_seed]
    if not U:
      return []

    ulo = min(l for l, r in U)

    # window families
    if not alt:
      # primary windows have a tiny deterministic shift per pass_id
      shift = (pass_id % 3) * 0.015
      window_fracs = [
        (0.12 + shift, 0.22 + shift),
        (0.35 + shift, 0.45 + shift),
        (0.58 + shift, 0.68 + shift),
        (0.80 + shift, 0.90 + shift),
      ]
    else:
      # alternate windows are offset to avoid stacking
      window_fracs = [
        (0.05, 0.15),
        (0.28, 0.38),
        (0.60, 0.70),
        (0.82, 0.92),
      ]

    # clamp windows inside (0.03, 0.97)
    window_fracs = [(max(0.03, a), min(0.97, b)) for (a, b) in window_fracs]

    blocks = []
    for idx, (fa, fb) in enumerate(window_fracs):
      win_lo = glo + int(round(fa * G))
      base = win_lo - ulo
      block = [(l + base, r + base) for (l, r) in U]
      # occasionally reverse internally to break symmetry deterministically
      if (idx + pass_id) % 2 == 1:
        block = list(reversed(block))
      blocks.append(block)

    # Interleave blocks to maximize mixing; pass_id determines order
    micro = []
    maxlen = max(len(b) for b in blocks)
    order = list(range(len(blocks)))
    if pass_id % 2 == 1:
      order.reverse()
    for i in range(maxlen):
      for j in order:
        blk = blocks[j]
        if i < len(blk):
          micro.append(blk[i])

    # micro-scale connectors (to tie windows)
    micro_connectors = [
      (glo + int(round(0.08 * G)), glo + int(round(0.30 * G))),
      (glo + int(round(0.60 * G)), glo + int(round(0.92 * G))),
      (glo + int(round(0.26 * G)), glo + int(round(0.56 * G))),
      (glo + int(round(0.44 * G)), glo + int(round(0.78 * G))),
    ]
    for a, b in micro_connectors:
      if b > a:
        micro.append((a, b))

    # alternate pass gets a long cross connector to tie distant colors
    if alt:
      long_a = glo + int(round(0.18 * G))
      long_b = glo + int(round(0.84 * G))
      if long_b > long_a:
        micro.append((long_a, long_b))

    # Trim to budget deterministically (take the prefix)
    if len(micro) > budget:
      micro = micro[:budget]
    return micro

  # Execute micro-phase A (primary)
  microA = build_micro_pass(T, passA_budget, pass_id=0, alt=False)
  if microA:
    take = min(len(microA), CAP - len(T))
    if take > 0:
      T.extend(microA[:take])

  # Optionally execute micro-phase B (alternate)
  if enable_alt_microphase:
    room = CAP - len(T)
    if room > 8:
      microB = build_micro_pass(T, passB_budget, pass_id=1, alt=True)
      if microB:
        take = min(len(microB), CAP - len(T))
        if take > 0:
          T.extend(microB[:take])

  # Final long-range deterministic connectors (4 small connectors)
  if len(T) < CAP:
    lo, hi, delta = span_info(T)
    span = max(1, delta)
    # Choose deterministic fractional anchors
    anchors = [0.02, 0.12, 0.28, 0.44]
    for a_frac in anchors:
      a = lo + int(round(a_frac * span))
      b = hi - int(round((a_frac + 0.01) * span))  # slightly asymmetric to avoid exact sym
      if b > a and len(T) < CAP:
        T.append((a, b))

  # Final tower-pins: short staggered intervals in interior windows
  room = CAP - len(T)
  if room > 8:
    lo, hi, delta = span_info(T)
    G = max(1, delta)
    eps = max(1, G // 4096)  # very short pins
    windows = [0.14, 0.31, 0.66, 0.83]
    # deterministic per-window count based on remaining capacity
    per_win = max(1, min(64, room // len(windows)))
    stride = max(2 * eps + 2, 3)
    pins = []
    for wi, w in enumerate(windows):
      base = lo + int(round(w * G))
      # start offset deterministically from wi to avoid alignment
      start = base - (per_win // 2) * stride + (wi % 3)
      for j in range(per_win):
        L = start + j * stride
        R = L + eps
        if R > L:
          pins.append((L, R))
    # append deterministic prefix up to room
    if pins:
      if len(pins) > room:
        pins = pins[:room]
      T.extend(pins)

  # Final trim to CAP
  if len(T) > CAP:
    T = T[:CAP]

  return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()