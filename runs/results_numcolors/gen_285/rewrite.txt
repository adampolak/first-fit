# EVOLVE-BLOCK-START

def construct_intervals(rounds=6,
                        rotate_starts=True,
                        reverse_block_parity=True,
                        interleave_blocks=True,
                        phase2_iters=1,
                        cross4_enabled=True):
  """
  Rotating-template KT spine with per-round fractional-window micro phases and
  final long-range connectors. Returns intervals (l, r) in FirstFit order.
  """

  # Capacity guard to keep total intervals < 10000
  CAP = 9800

  # Deterministic seed for reproducibility of choices
  BASE_SEED = 1729

  # Start templates (Kierstead–Trotter–style)
  spine_starts = (2, 6, 10, 14)
  template_bank = [
    (2, 6, 10, 14),  # classic KT
    (1, 5, 9, 13),   # left-shifted
    (3, 7, 11, 15),  # right-shifted
    (4, 8, 12, 16),  # stretched-right
  ]

  # Seed: a single short interval
  T = [(0, 1)]

  # Helpers
  def _span(current_T):
    lo = min(l for l, r in current_T)
    hi = max(r for l, r in current_T)
    delta = hi - lo
    if delta <= 0:
      delta = 1
    return lo, hi, delta

  def _intv(l, r):
    li, ri = int(l), int(r)
    if ri <= li:
      ri = li + 1
    return (li, ri)

  def _seed_pick(seed, ridx, salt=0):
    # Simple deterministic mixing
    return (hash((BASE_SEED, seed, ridx, salt)) & 0x7fffffff)

  def _build_blocks(base_T, starts, lo, delta, reverse_parity=False):
    """
    Build translated blocks; optionally reverse every odd block to mix FF colors.
    """
    blocks = []
    for b_idx, s in enumerate(starts):
      block_src = base_T[::-1] if (reverse_parity and (b_idx % 2 == 1)) else base_T
      base = s * delta - lo
      block = [(l + base, r + base) for (l, r) in block_src]
      blocks.append(block)
    return blocks

  def _interleave_blocks(blocks, ridx, interleave=True):
    """
    Round-robin interleaving to maximize mixing; parity flip by ridx.
    """
    if not interleave:
      out = []
      for blk in blocks:
        out.extend(blk)
      return out

    order = list(range(len(blocks)))
    if (ridx % 2) == 1:
      order = order[::-1]
    # small rotation to avoid identical prefix each round
    if order:
      krot = ridx % len(order)
      order = order[krot:] + order[:krot]

    maxlen = max((len(b) for b in blocks), default=0)
    out = []
    for i in range(maxlen):
      for idx in order:
        blk = blocks[idx]
        if i < len(blk):
          out.append(blk[i])
    return out

  def _append_classic_connectors(S, starts, delta):
    s0, s1, s2, s3 = starts
    S.append(((s0 - 1) * delta, (s1 - 1) * delta))  # left cap
    S.append(((s2 + 2) * delta, (s3 + 2) * delta))  # right cap
    S.append(((s0 + 2) * delta, (s2 - 1) * delta))  # cross 1
    S.append(((s1 + 2) * delta, (s3 - 1) * delta))  # cross 2

  def _per_round_micro_layer(pre_round_T, ridx, whole_round_S, cap_room):
    """
    A small deterministic micro-phase after each round to densify the backbone
    without inflating omega. Uses a sparse seed and short pinned copies.
    """
    if cap_room <= 12 or not pre_round_T:
      return []

    glo = min(l for l, r in whole_round_S)
    ghi = max(r for l, r in whole_round_S)
    G = max(1, ghi - glo)

    # Thin, evenly spaced seed from pre-round T
    seed_sz = max(10, min(28, len(pre_round_T) // 350 if len(pre_round_T) >= 350 else 14))
    stride = max(1, len(pre_round_T) // max(1, seed_sz))
    U = [pre_round_T[i] for i in range(0, len(pre_round_T), stride)][:seed_sz]
    if not U:
      return []

    ulo = min(l for l, r in U)

    # Two window families; pick by deterministic seed from BASE_SEED and ridx
    windows_A = [(0.06, 0.14), (0.28, 0.38), (0.54, 0.64), (0.78, 0.88)]
    windows_B = [(0.12, 0.22), (0.35, 0.45), (0.58, 0.68), (0.80, 0.90)]
    pickB = (_seed_pick(7, ridx) % 3 == 1)
    windows = windows_B if pickB else windows_A

    # Use short pins to limit clique growth: pin length ~ G/256 (at least 1)
    pin_len = max(1, G // 256)

    micro_blocks = []
    for (fa, fb) in windows:
      win_lo = glo + int(round(fa * G))
      base = win_lo - ulo
      # Pin each translated copy near its midpoint
      block = []
      for idx, (l, r) in enumerate(U):
        mid = (l + r) // 2
        L = mid + base - (pin_len // 2) + (idx % 2)  # tiny stagger
        R = L + pin_len
        if R > L:
          block.append((L, R))
      # Reverse every other micro-block to further break symmetry
      if len(micro_blocks) % 2 == 1:
        block = list(reversed(block))
      micro_blocks.append(block)

    # Interleave micro-blocks (reverse order on odd rounds)
    micro = []
    maxlen = max((len(b) for b in micro_blocks), default=0)
    reverse_order = ((ridx % 2) == 1)
    for i in range(maxlen):
      for blk in (reversed(micro_blocks) if reverse_order else micro_blocks):
        if i < len(blk):
          micro.append(blk[i])

    # Fractional-span micro connectors (short) to tie the windows
    micro_connectors = [
      (glo + int(round(0.08 * G)), glo + int(round(0.08 * G)) + pin_len),
      (glo + int(round(0.32 * G)), glo + int(round(0.32 * G)) + pin_len),
      (glo + int(round(0.66 * G)), glo + int(round(0.66 * G)) + pin_len),
      (glo + int(round(0.90 * G)), glo + int(round(0.90 * G)) + pin_len),
    ]
    for a, b in micro_connectors:
      if b > a:
        micro.append((a, b))

    # Trim to available room (keep conservative per round)
    take = min(len(micro), max(0, min(cap_room, 64)))
    return micro[:take]

  def _final_long_connectors(T, cap_room, enable=True):
    """
    Add a very small set of long-range connectors spanning 4–6 deltas
    to tie distant scales. Conservative to avoid omega blow-up.
    """
    if not enable or cap_room <= 4 or not T:
      return []
    lo = min(l for l, r in T)
    hi = max(r for l, r in T)
    G = max(1, hi - lo)
    # Two long cross connectors and two side caps
    L1 = (lo + int(round(0.14 * G)), lo + int(round(0.86 * G)))
    L2 = (lo + int(round(0.22 * G)), lo + int(round(0.78 * G)))
    C1 = (lo + int(round(0.05 * G)), lo + int(round(0.35 * G)))
    C2 = (lo + int(round(0.65 * G)), lo + int(round(0.95 * G)))
    longs = [L1, L2, C1, C2]
    out = []
    for a, b in longs:
      if b > a:
        out.append((a, b))
      if len(out) >= cap_room:
        break
    return out

  # Stage 1: KT spine with parity interleaving and per-round micro-layers
  rounds = max(0, int(rounds))
  for ridx in range(rounds):
    # Size growth prediction: sz -> 4*sz + 4 (classic connectors only),
    # micro-layer adds at most ~64 intervals per round (CAP-guarded).
    if 4 * len(T) + 4 > CAP:
      break

    # Snapshot pre-round spine for micro-layer construction
    pre_T = T

    # Choose starts
    starts = template_bank[ridx % len(template_bank)] if rotate_starts else spine_starts

    # Compute geometry from pre-round T
    lo, hi, delta = _span(pre_T)

    # Build four translated blocks with optional reverse parity
    reverse_parity = bool(reverse_block_parity and ((ridx % 2) == 1))
    blocks = _build_blocks(pre_T, starts, lo, delta, reverse_parity=reverse_parity)

    # Interleave blocks using deterministic parity schedule
    do_inter = bool(interleave_blocks and ((ridx % 2) == 0))
    S = _interleave_blocks(blocks, ridx, interleave=do_inter)

    # Append classic KT connectors
    _append_classic_connectors(S, starts, delta)

    # Capacity-aware per-round micro-layer densification
    room_now = CAP - len(S)
    micro = _per_round_micro_layer(pre_T, ridx, S, room_now)
    if micro:
      S.extend(micro)

    # Update T after round
    T = S

    # If we are very close to capacity, stop early
    if len(T) >= CAP - 64:
      break

  # Global micro-phase(s) – retain compatibility with phase2_iters
  def build_micro_delta_round(current_T, budget, ridx_seed):
    if not current_T or budget <= 8:
      return []
    glo = min(l for l, r in current_T)
    ghi = max(r for l, r in current_T)
    G = max(1, ghi - glo)

    # Thin seed
    seed_sz = max(12, min(32, len(current_T) // 300 if len(current_T) >= 300 else 24))
    stride = max(1, len(current_T) // max(1, seed_sz))
    U = [current_T[i] for i in range(0, len(current_T), stride)][:seed_sz]
    if not U:
      return []
    ulo = min(l for l, r in U)

    # Alternate window families deterministically from ridx_seed
    windows_A = [(0.08, 0.18), (0.30, 0.40), (0.58, 0.68), (0.80, 0.90)]
    windows_B = [(0.05, 0.15), (0.28, 0.38), (0.60, 0.70), (0.82, 0.92)]
    use_B = (_seed_pick(11, ridx_seed) % 2 == 0)
    windows = windows_B if use_B else windows_A

    blocks = []
    for (fa, fb) in windows:
      win_lo = glo + int(round(fa * G))
      base = win_lo - ulo
      block = [(l + base, r + base) for (l, r) in U]
      blocks.append(block)

    # Interleave with reversed order to diversify coupling
    micro = []
    maxlen = max((len(b) for b in blocks), default=0)
    for i in range(maxlen):
      for blk in reversed(blocks):
        if i < len(blk):
          micro.append(blk[i])

    # Micro connectors (medium length)
    micro_connectors = [
      (glo + int(round(0.10 * G)), glo + int(round(0.32 * G))),
      (glo + int(round(0.48 * G)), glo + int(round(0.72 * G))),
      (glo + int(round(0.22 * G)), glo + int(round(0.56 * G))),
      (glo + int(round(0.66 * G)), glo + int(round(0.90 * G))),
    ]
    for a, b in micro_connectors:
      if b > a:
        micro.append((a, b))

    if len(micro) > budget:
      micro = micro[:budget]
    return micro

  steps = min(max(0, int(phase2_iters)), 2)
  for step in range(steps):
    room = CAP - len(T)
    if room <= 8:
      break
    micro_global = build_micro_delta_round(T, room, ridx_seed=(rounds + step))
    if micro_global:
      if len(micro_global) > room:
        micro_global = micro_global[:room]
      T.extend(micro_global)

  # Final long-range deterministic connectors (cross-scale)
  room = CAP - len(T)
  if room > 0:
    longs = _final_long_connectors(T, room=min(room, 4), enable=bool(cross4_enabled))
    if longs:
      T.extend(longs)

  # Normalize to non-negative integer endpoints with positive length and cap to CAP
  if not T:
    return []
  min_l = min(l for l, r in T)
  if min_l < 0:
    T = [(l - min_l, r - min_l) for (l, r) in T]
  out = [_intv(l, r) for (l, r) in T]
  if len(out) > CAP:
    out = out[:CAP]
  return out

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()