# EVOLVE-BLOCK-START

def construct_intervals(seed_count=1, rounds=6, BASE_SEED=123457):
  """
  Deterministic spine + layered micro-phases interval construction.

  Signature preserved: construct_intervals(seed_count=1) -> list of (l, r) tuples.
  """

  CAP = 9800  # hard capacity guard (< 10000)
  # Template bank for four-block spine translations (rotating)
  template_bank = [
    (2, 6, 10, 14),
    (1, 5, 9, 13),
    (3, 7, 11, 15),
    (4, 8, 12, 16),
  ]

  # Seed intervals (keep small to allow many rounds)
  if seed_count <= 1:
    T = [(0, 1)]
  else:
    step = 3
    T = [(i * step, i * step + 1) for i in range(min(seed_count, 8))]

  # --- Utility helpers (modularized) ---
  def _span(current):
    lo = min(l for l, _ in current)
    hi = max(r for _, r in current)
    delta = hi - lo
    if delta <= 0:
      delta = 1
    return lo, hi, delta

  def _translate_block(current, base):
    return [(l + base, r + base) for (l, r) in current]

  def _interleave_blocks(blocks, order=None):
    if order is None:
      order = list(range(len(blocks)))
    maxlen = max((len(b) for b in blocks), default=0)
    out = []
    for i in range(maxlen):
      for idx in order:
        blk = blocks[idx]
        if i < len(blk):
          out.append(blk[i])
    return out

  def _append_connectors(seq, starts, delta, lo, add_long=False):
    # classic four connectors (keep endpoints inside the coarse grid)
    s0, s1, s2, s3 = starts
    # compute connector endpoints relative to lo (so we stay consistent with block translations)
    seq.append(((s0 - 1) * delta + lo, (s1 - 1) * delta + lo))  # left cap
    seq.append(((s2 + 2) * delta + lo, (s3 + 2) * delta + lo))  # right cap
    seq.append(((s0 + 2) * delta + lo, (s2 - 1) * delta + lo))  # cross 1
    seq.append(((s1 + 2) * delta + lo, (s3 - 1) * delta + lo))  # cross 2
    if add_long:
      # one cautious long cross to tie far ends (span ~4*delta)
      seq.append(((s0 + 4) * delta + lo, (s3 + 4) * delta + lo))

  def _cap_inside(lo, hi, a_frac, b_frac):
    span = max(1, hi - lo)
    L = lo + max(1, int(round(a_frac * span)))
    R = lo + max(1, int(round(b_frac * span)))
    if R <= L:
      R = L + 1
    return (L, R)

  def _duplicate_with_offset(blocks, delta, lo, offset, cap_left):
    """
    Produce an extra copy of each block shifted by 'offset'.
    Stop early if adding all would exceed cap_left.
    Returns list of extra intervals (in order).
    """
    extras = []
    for blk in blocks:
      for (l, r) in blk:
        if cap_left <= 0:
          return extras
        nl, nr = l + offset, r + offset
        # keep interval valid
        if nr > nl:
          extras.append((nl, nr))
          cap_left -= 1
    return extras

  # deterministic per-round "random" choices derived from BASE_SEED
  def _per_round_seed(round_idx):
    # simple deterministic mixing
    return (BASE_SEED * 1000003 + round_idx * 1009) & 0xFFFFFFFF

  # --- Spine construction (modular KT-style rounds with CAP-aware density) ---
  for ridx in range(max(0, int(rounds))):
    # predictive guard: size -> 4*size + 4 (+ possible duplicates)
    predicted = 4 * len(T) + 4
    if predicted > CAP:
      break

    starts = template_bank[ridx % len(template_bank)]
    lo, hi, delta = _span(T)
    # compute block bases
    blocks = []
    for s in starts:
      base = s * delta + lo - lo  # keep base as s*delta (consistent with connector math)
      # Since we translate by (s*delta - lo), make base = s*delta - lo + lo = s*delta
      # (we will add lo later when computing connectors). To keep integers, set base = s*delta - lo_offset
      base = s * delta - lo
      blocks.append(_translate_block(T, base))

    # deterministic interleave flag toggles every round (parity)
    per_seed = _per_round_seed(ridx)
    do_interleave = (ridx % 2 == 0)  # parity schedule (reproducible)
    # order reversal for parity mixing
    order = list(range(4))
    if not do_interleave:
      # sequential order, sometimes reversed to break symmetry
      if (per_seed >> 3) & 1:
        order = list(reversed(order))
      out_seq = []
      for idx in order:
        out_seq.extend(blocks[idx])
    else:
      # interleave in deterministic order (rotate order by low bits)
      rot = (per_seed & 3) % 4
      inter_order = [(i + rot) % 4 for i in range(4)]
      out_seq = _interleave_blocks(blocks, order=inter_order)

    # CAP-aware density multiplier: on even rounds attempt one duplicate copy with a small offset
    density_allowed = (len(T) + len(out_seq) + 4) < CAP  # basic allowance before extras
    extras = []
    if density_allowed:
      # choose small integer offset deterministic from seed (1..max_off)
      max_off = max(1, delta // 8)
      chosen_off = 1 + ((per_seed >> 5) % max_off)
      # but only attempt duplication if we would not exceed CAP; duplicate only once
      cap_left = CAP - (len(T) + len(out_seq) + 4)
      extras = _duplicate_with_offset([out_seq], delta, lo, chosen_off, cap_left)
      # extras is a list; note we used out_seq packaged as single block list above

    # append backbone round blocks to T (order preserved)
    T = out_seq.copy()
    # append connectors for this round (use add_long only on final executed spine round)
    add_long = (ridx == max(0, int(rounds)) - 1)
    _append_connectors(T, starts, delta, lo, add_long=add_long)

    # append extras (density duplicates) after connectors to push FirstFit mixing
    if extras:
      # extras were generated as copies of out_seq; append until cap
      room = CAP - len(T)
      if room > 0:
        if len(extras) > room:
          extras = extras[:room]
        T.extend(extras)

  # If near capacity, return early
  if len(T) >= CAP - 24:
    return T[:CAP]

  # --- Layered micro-phases (two distinct window families) ---
  # Build a thin seed U from T (evenly sampled) for micro-block translations
  def _thin_seed(current, target_sz):
    n = len(current)
    if n == 0:
      return []
    step = max(1, n // max(1, target_sz))
    return [current[i] for i in range(0, n, step)][:target_sz]

  lo_g, hi_g, delta_g = _span(T)

  # two micro window families
  window_fracs1 = [(0.12, 0.22), (0.35, 0.45), (0.58, 0.68), (0.80, 0.90)]
  window_fracs2 = [(0.06, 0.14), (0.28, 0.38), (0.54, 0.64), (0.78, 0.88)]

  # micro-builder shared logic
  def build_micro(current_T, window_fracs, round_tag=0, alt_blockers=False):
    room = CAP - len(current_T)
    if room <= 8 or not current_T:
      return []

    glo = min(l for l, r in current_T)
    ghi = max(r for l, r in current_T)
    G = max(1, ghi - glo)

    # choose thin seed size deterministically (bounded)
    target_seed =  max(8, min(48, len(current_T) // 250))
    U = _thin_seed(current_T, target_seed)
    if not U:
      return []

    ulo = min(l for l, r in U)

    # build blocks translated into each window
    blocks = []
    for (fa, fb) in window_fracs:
      win_lo = glo + int(round(fa * G))
      base = win_lo - ulo
      block = [(l + base, r + base) for (l, r) in U]
      blocks.append(block)

    # optionally create small blocker pins before blocks (alternate pass)
    pins = []
    if alt_blockers:
      # create short pins inside each window to occupy low-index colors early
      for (fa, fb) in window_fracs:
        w_lo = glo + int(round(fa * G))
        w_hi = glo + int(round(fb * G))
        width = max(2, w_hi - w_lo)
        pin_len = max(1, max(1, width // 12))
        # create a handful of staggered pins
        for i in range(3):
          L = w_lo + (i * max(1, (width - pin_len))) // 2
          R = min(w_hi, L + pin_len)
          if R > L:
            pins.append((L, R))

    # interleave blocks (order depends on round_tag for diversity)
    order = list(range(len(blocks)))
    if (round_tag % 2) == 1:
      order.reverse()
    micro_body = _interleave_blocks(blocks, order=order)

    # add small micro-connectors that tie the windows together
    micro_conn = [
      (glo + int(round(0.08 * G)), glo + int(round(0.30 * G))),
      (glo + int(round(0.60 * G)), glo + int(round(0.92 * G))),
      (glo + int(round(0.26 * G)), glo + int(round(0.56 * G))),
    ]

    out = []
    # put pins first if any (to occupy small colors early)
    if pins:
      out.extend(pins)
    out.extend(micro_body)
    out.extend([c for c in micro_conn if c[1] > c[0]])

    # trim to available room
    if len(out) > room:
      out = out[:room]
    return out

  # Primary micro-phase using window_fracs1, alt_blockers False
  room = CAP - len(T)
  if room > 16:
    micro1 = build_micro(T, window_fracs1, round_tag=0, alt_blockers=False)
    if micro1:
      room_now = CAP - len(T)
      if len(micro1) > room_now:
        micro1 = micro1[:room_now]
      T.extend(micro1)

  # Secondary micro-phase using window_fracs2 with blockers
  room = CAP - len(T)
  if room > 12:
    micro2 = build_micro(T, window_fracs2, round_tag=1, alt_blockers=True)
    if micro2:
      room_now = CAP - len(T)
      if len(micro2) > room_now:
        micro2 = micro2[:room_now]
      T.extend(micro2)

  # Cross-scale connectors: deterministic long-range connectors spanning ~4-6 delta_g
  room = CAP - len(T)
  if room > 0 and T:
    lo, hi, delta = _span(T)
    # anchor connectors to spine geometry: pick few offsets in {4,5,6} * delta
    long_spans = [4, 5, 6]
    connectors = []
    # compute anchors deterministically
    per_seed = _per_round_seed(9999)
    for i, ls in enumerate(long_spans):
      if len(connectors) >= 6:
        break
      # choose anchors using low bits
      a_offset = (i * 3 + ((per_seed >> (i + 2)) & 3)) % 8
      b_offset = a_offset + ls
      L = lo + a_offset * delta // max(1, 2)
      R = lo + b_offset * delta // max(1, 2)
      if R > L:
        connectors.append((L, R))
    # append until capacity
    for c in connectors:
      if CAP - len(T) <= 0:
        break
      T.append(c)

  # Final lightweight cap-trio in parabolic placement (controlled, capacity-guarded)
  room = CAP - len(T)
  if room > 3:
    lo, hi, delta = _span(T)
    cap1 = _cap_inside(lo, hi, 0.15, 0.60)
    cap2 = _cap_inside(lo, hi, 0.25, 0.80)
    cap3 = _cap_inside(lo, hi, 0.55, 0.95)
    for c in (cap1, cap2, cap3):
      if CAP - len(T) <= 0:
        break
      # be conservative: only add if it extends the range meaningfully
      if c[1] - c[0] >= max(1, delta // 20):
        T.append(c)

  # Final trim to CAP and return
  if len(T) > CAP:
    T = T[:CAP]
  return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()