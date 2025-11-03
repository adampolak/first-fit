# EVOLVE-BLOCK-START

def _span_delta(T):
  lo = min(l for l, _ in T)
  hi = max(r for _, r in T)
  delta = hi - lo
  if delta <= 0:
    delta = 1
  return lo, hi, delta

def _per_round_seed(base_seed, idx):
  # Deterministic per-round integer derived from base seed and index
  # Small LCG-like mix for reproducibility without importing random/hash
  x = (base_seed ^ (idx * 0x9e3779b1)) & 0xffffffff
  x = (x * 1664525 + 1013904223) & 0xffffffff
  return x

def _build_blocks(current_T, starts, lo, delta, do_interleave=False, reverse_order=False, K=1):
  """
  Build translated blocks from current_T using spacing multiplier K.
  Interleaving and block reversal options supported.
  """
  blocks = []
  for s in starts:
    base = s * delta * K - lo
    block = [(l + base, r + base) for (l, r) in current_T]
    blocks.append(block)

  # Assemble result according to interleaving flag
  S = []
  if do_interleave:
    order = list(range(len(blocks)))
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
  return S

def _append_classic_connectors(S, starts, delta, add_cross4=False, long_offset=None):
  """
  Append classic four connectors (caps + crosses). Optionally add a long-range cross4
  and an explicit long_offset (span multiplier) connectors list.
  """
  s0, s1, s2, s3 = starts
  S.append(((s0 - 1) * delta, (s1 - 1) * delta))  # left cap
  S.append(((s2 + 2) * delta, (s3 + 2) * delta))  # right cap
  S.append(((s0 + 2) * delta, (s2 - 1) * delta))  # cross 1
  S.append(((s1 + 2) * delta, (s3 - 1) * delta))  # cross 2
  if add_cross4:
    S.append(((s0 + 4) * delta, (s3 + 4) * delta))
  if long_offset is not None:
    # anchor at s0 and s3 span long_offset * delta
    S.append(((s0 + long_offset) * delta, (s3 + long_offset) * delta))

def _insert_near_tail(seq, intervals):
  """
  Insert intervals near the tail staggered so each new interval intersects many
  active colors late in the presentation.
  """
  out = list(seq)
  for i, iv in enumerate(intervals):
    pos = len(out) - (i * 2 + 1)
    if pos < 0:
      out.append(iv)
    else:
      out.insert(pos, iv)
  return out

def build_micro_round(current_T, budget, iter_id, base_seed, window_fracs):
  """
  Build a delta-scale micro round from a thin seed of current_T.
  - budget: maximum intervals to produce (trimmed)
  - window_fracs: list of (fa, fb) fractions to place blocks
  - returns list of intervals (already capacity-trimmed)
  """
  if not current_T or budget <= 6:
    return []

  glo = min(l for l, r in current_T)
  ghi = max(r for l, r in current_T)
  G = max(1, ghi - glo)

  # Thin evenly-spaced seed (bounded)
  max_seed = 36
  stride = max(1, len(current_T) // max_seed)
  U = [current_T[i] for i in range(0, len(current_T), stride)][:max_seed]
  if not U:
    return []

  ulo = min(l for l, r in U)

  # Deterministic tag to vary minor choices
  tag = _per_round_seed(base_seed, iter_id)
  do_interleave = (tag % 2 == 0)  # alternate interleaving per micro round deterministically
  reverse_block_order = ((tag >> 1) % 2 == 1)

  blocks = []
  for bi, (fa, fb) in enumerate(window_fracs):
    win_lo = glo + int(round(max(0.01, min(0.98, fa)) * G))
    base = win_lo - ulo
    block = [(l + base, r + base) for (l, r) in U]
    # Alternate internal reversal to break symmetry
    if ((bi + tag) % 2) == 1:
      block = list(reversed(block))
    blocks.append(block)

  # Interleave micro-blocks to maximize mixing
  packed = []
  maxlen = max(len(b) for b in blocks)
  order = list(range(len(blocks)))
  if reverse_block_order:
    order.reverse()
  for i in range(maxlen):
    for idx in order:
      blk = blocks[idx]
      if i < len(blk):
        packed.append(blk[i])

  # Add micro connectors (fractional caps and crosses)
  connectors = []
  connectors.append((glo + int(round(0.08 * G)), glo + int(round(0.30 * G))))
  connectors.append((glo + int(round(0.60 * G)), glo + int(round(0.92 * G))))
  connectors.append((glo + int(round(0.26 * G)), glo + int(round(0.56 * G))))
  connectors.append((glo + int(round(0.44 * G)), glo + int(round(0.78 * G))))

  # Optionally add a longer micro cross anchored across multiple windows (choose 4-6 delta-span analog)
  long_span_choice = 4 + (tag % 3)  # 4,5,6 cyclic
  # Convert to micro-scale by mapping offsets into global coordinates:
  # anchor near window extremes to cross many windows
  a = glo + int(round(0.05 * G))
  b = glo + int(round(0.95 * G))
  connectors.append((a, b))  # big micro connector

  # Compose packed followed by connectors
  out = packed + [c for c in connectors if c[1] > c[0]]

  # Trim to budget
  if len(out) > budget:
    out = out[:budget]
  return out

def construct_intervals(seed_count=1):
  """
  Improved KT-spine with deterministic parity interleaving, CAP-aware density multiplier,
  and two distinct micro-phases with long-range cross connectors. Returns list of (l, r)
  in FirstFit presentation order. Maintains the same interface as prior versions.
  """
  CAP = 9800
  BASE_SEED = 1234567  # deterministic seed for per-round choices

  template_bank = [
    (2, 6, 10, 14),  # T1 classic KT
    (1, 5, 9, 13),   # T2 left-shifted
    (3, 7, 11, 15),  # T3 right-shifted
    (4, 8, 12, 16),  # T4 stretched right
  ]

  # Seed with a small deterministic multi-seed if requested, otherwise single unit
  if seed_count <= 1:
    T = [(0, 1)]
  else:
    # deterministic sparse seeds up to 4
    sc = max(1, int(seed_count))
    sc = min(sc, 4)
    step = 2
    T = [(i * step, i * step + 1) for i in range(sc)]

  # Main spine: up to 6 rounds but guard by CAP. Use alternating parity: even rounds interleave.
  MAX_ROUNDS = 6
  for ridx in range(MAX_ROUNDS):
    # Predict conservative next size and stop early if would blow CAP
    predicted = 4 * len(T) + 4 + 8  # +8 to account for typical connectors
    if predicted > CAP:
      break

    lo, hi, delta = _span_delta(T)
    starts = template_bank[ridx % len(template_bank)]

    # Deterministic per-round choices
    tag = _per_round_seed(BASE_SEED, ridx)
    do_interleave = (ridx % 2 == 0)  # parity schedule: even rounds interleave (deterministic)
    reverse_order = (ridx % 2 == 1)  # odd rounds reverse block order to break repetition

    # CAP-aware small density multiplier K: occasionally increase spacing to 2 to couple colors
    # Use per-round tag to decide limited K=2 on some rounds (but ensure budget)
    K = 2 if (tag % 5 == 0 and len(T) > 16) else 1  # deterministic sparsely applied

    S = _build_blocks(T, starts, lo, delta, do_interleave=do_interleave, reverse_order=reverse_order, K=K)

    # Append classic connectors; add an extra cross4 on final executed round if budget allows
    add_cross4 = (ridx == MAX_ROUNDS - 1)
    # Choose a long_offset (4..6) deterministically to add a long-range connector sometimes
    long_offset = 4 + (tag % 3) if (tag % 7 == 0) else None
    _append_classic_connectors(S, starts, delta, add_cross4=add_cross4, long_offset=long_offset)

    # If adding S would exceed CAP, trim S conservatively (keep head since earlier items are more influential)
    room = CAP - len(S)
    if room < 0:
      # current S alone exceeds cap (rare): truncate S
      S = S[:CAP]
      T = S
      break

    # Replace T by S (spine growth)
    T = S
    if len(T) >= CAP:
      T = T[:CAP]
      break

  # If already near capacity, return
  if len(T) >= CAP - 16:
    return T

  # Micro-phase 1: fractional windows (primary)
  room = CAP - len(T)
  if room > 8:
    window_fracs = [(0.12, 0.22), (0.35, 0.45), (0.58, 0.68), (0.80, 0.90)]
    micro1 = build_micro_round(T, room // 2, iter_id=0, base_seed=BASE_SEED + 11, window_fracs=window_fracs)
    # Insert micro1 near the tail to intersect many colors
    if micro1:
      avail = CAP - len(T)
      micro1 = micro1[:avail]
      T = _insert_near_tail(T, micro1)

  # Micro-phase 2: alternate windows (secondary)
  room = CAP - len(T)
  if room > 8:
    window_fracs2 = [(0.06, 0.14), (0.28, 0.38), (0.54, 0.64), (0.78, 0.88)]
    micro2 = build_micro_round(T, room, iter_id=1, base_seed=BASE_SEED + 17, window_fracs=window_fracs2)
    if micro2:
      avail = CAP - len(T)
      micro2 = micro2[:avail]
      T = _insert_near_tail(T, micro2)

  # Final small tail: three long symmetric caps inserted in tail positions to force late FF colors
  room = CAP - len(T)
  if room > 0:
    lo, hi, delta = _span_delta(T)
    if delta <= 0:
      delta = 1
    # Build caps but ensure viable endpoints
    caps = [
      (lo + max(1, int(round(0.08 * delta))), lo + max(1, int(round(0.60 * delta)))),
      (lo + max(1, int(round(0.25 * delta))), lo + max(1, int(round(0.75 * delta)))),
      (lo + max(1, int(round(0.45 * delta))), lo + max(1, int(round(0.92 * delta))))
    ]
    caps = [c for c in caps if c[1] > c[0]]
    # Insert only up to room intervals
    to_add = caps[:room]
    if to_add:
      T = _insert_near_tail(T, to_add)

  # Final safety trim
  if len(T) > CAP:
    T = T[:CAP]
  return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()