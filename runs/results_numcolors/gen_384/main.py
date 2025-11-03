# EVOLVE-BLOCK-START

def _span_delta(T):
  lo = min(l for l, _ in T)
  hi = max(r for _, r in T)
  delta = hi - lo
  if delta <= 0:
    delta = 1
  return lo, hi, delta

def _append_connectors(S, delta, starts):
  # Classic four connectors providing cross-links between translated blocks
  s0, s1, s2, s3 = starts
  S.append(((s0 - 1) * delta, (s1 - 1) * delta))  # left cap
  S.append(((s2 + 2) * delta, (s3 + 2) * delta))  # right cap
  S.append(((s0 + 2) * delta, (s2 - 1) * delta))  # cross 1
  S.append(((s1 + 2) * delta, (s3 - 1) * delta))  # cross 2

def _apply_round(current_T, starts, do_interleave=False, reverse_order=False):
  lo, hi, delta = _span_delta(current_T)

  # Build four translated blocks; deterministic, compact construction
  blocks = []
  for s in starts:
    base = s * delta - lo
    block = [(l + base, r + base) for (l, r) in current_T]
    blocks.append(block)

  # Assemble S with optional interleaving
  S = []
  if do_interleave:
    maxlen = max(len(b) for b in blocks)
    order = [0, 1, 2, 3]
    if reverse_order:
      order = list(reversed(order))
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

  _append_connectors(S, delta, starts)
  return S

def construct_intervals(rounds=6,
                        rotate_starts=True,
                        reverse_block_parity=True,
                        interleave_blocks=True,
                        phase2_iters=2):
  """
  Produces a deterministic sequence of intervals with a rotating KT-spine backbone
  and a small micro-phase, designed to fit within CAP=9800.
  Returns: list[(l, r)]
  """

  CAP = 9800
  spine_starts = [2, 6, 10, 14]

  # Rotating template bank (kept small for speed and stability)
  template_bank = [
    [2, 6, 10, 14],  # T1 classic
    [1, 5, 9, 13],   # T2 left-shifted
    [3, 7, 11, 15],  # T3 right-shifted
    [4, 8, 12, 16],  # T4 stretched-right
  ]

  # Seed with a single unit interval
  T = [(0, 1)]

  # Helper: cap-aware round growth
  def next_round_size(sz):
    return 4 * sz + 4

  # Determine feasible depth (depth <= rounds) within CAP
  depth = 0
  cur = len(T)
  while depth < rounds:
    nxt = next_round_size(cur)
    if nxt > CAP:
      break
    cur = nxt
    depth += 1

  # Stage 1: six KT rounds; interleave/parity controlled by flags
  for ridx in range(depth):
    starts = template_bank[ridx % len(template_bank)] if rotate_starts else spine_starts
    do_interleave = bool(interleave_blocks and (ridx % 2 == 0))
    rev = bool(reverse_block_parity and (ridx % 2 == 1))
    T = _apply_round(T, starts, do_interleave=do_interleave, reverse_order=rev)
    if len(T) >= CAP:
      return T[:CAP]

  if len(T) >= CAP - 16:
    return T[:CAP]

  # Stage 2: compact micro-phase
  def cap_at(a_frac, b_frac, lo, hi):
    delta = hi - lo
    if delta <= 0:
      delta = 1
    L = lo + max(1, int(round(a_frac * delta)))
    R = lo + max(1, int(round(b_frac * delta)))
    if R <= L:
      R = L + 1
    return (L, R)

  lo, hi, delta = _span_delta(T)

  tail_caps = [
    cap_at(0.08, 0.60, lo, hi),
    cap_at(0.25, 0.75, lo, hi),
    cap_at(0.75, 0.92, lo, hi),
  ]

  # Insert tail caps in a way that preserves determinism and budget
  room = CAP - len(T)
  if room > 0:
    # deterministic minimal insertion from the tail
    for c in tail_caps:
      if len(T) >= CAP:
        break
      T.append(c)

  if len(T) >= CAP - 16:
    return T[:CAP]

  # Micro-phase B: delta-driven micro rounds
  def thin_seed(current_T, max_seed):
    n = len(current_T)
    if n == 0 or max_seed <= 0:
      return []
    step = max(1, n // max_seed)
    return current_T[::step][:max_seed]

  def build_micro_delta_round(current_T, budget, iter_id=0, alt=False):
    if not current_T or budget <= 8:
      return []

    glo = min(l for l, _ in current_T)
    ghi = max(r for _, r in current_T)
    G = max(1, ghi - glo)

    # Seed from current_T
    seed_sz = max(8, min(40, len(current_T) // 250))
    stride = max(1, len(current_T) // seed_sz)
    U = [current_T[i] for i in range(0, len(current_T), stride)][:seed_sz]
    if not U:
      return []

    ulo = min(l for l, _ in U)

    window_fracs = [(0.12, 0.22), (0.35, 0.45), (0.58, 0.68), (0.80, 0.90)]
    blocks = []
    for (fa, fb) in window_fracs:
      win_lo = glo + int(round(fa * G))
      base = win_lo - ulo
      block = [(l + base, r + base) for (l, r) in U]
      blocks.append(block)

    micro = []
    maxlen = max(len(b) for b in blocks)
    for i in range(maxlen):
      for blk in blocks:
        if i < len(blk):
          micro.append(blk[i])

    connectors = [
      (glo + int(round(0.08 * G)), glo + int(round(0.30 * G))),
      (glo + int(round(0.60 * G)), glo + int(round(0.92 * G))),
      (glo + int(round(0.26 * G)), glo + int(round(0.56 * G))),
      (glo + int(round(0.44 * G)), glo + int(round(0.78 * G))),
    ]
    if alt:
      connectors.append((glo + int(round(0.18 * G)), glo + int(round(0.84 * G))))
    for a, b in connectors:
      if b > a:
        micro.append((a, b))

    if len(micro) > budget:
      micro = micro[:budget]
    return micro

  # First micro pass
  remaining = CAP - len(T)
  if remaining > 8:
    micro1 = build_micro_delta_round(T, remaining, iter_id=0, alt=False)
    if micro1:
      if len(micro1) > remaining:
        micro1 = micro1[:remaining]
      T.extend(micro1)

  # Second micro pass
  remaining = CAP - len(T)
  if remaining > 0:
    micro2 = build_micro_delta_round(T, remaining, iter_id=1, alt=False)
    if micro2:
      if len(micro2) > remaining:
        micro2 = micro2[:remaining]
      T.extend(micro2)

  # Return with normalization to integers and non-degenerate intervals
  if not T:
    return []
  min_l = min(l for l, _ in T)
  if min_l < 0:
    shift = -min_l
    T = [(l + shift, r + shift) for (l, r) in T]

  intervals = []
  for (l, r) in T:
    li = int(round(l))
    ri = int(round(r))
    if ri <= li:
      ri = li + 1
    intervals.append((li, ri))

  return intervals


def run_experiment(**kwargs):
  """Main called by evaluator"""
  # Forward all kwargs to preserve compatibility with the original interface
  return construct_intervals(**kwargs)

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()