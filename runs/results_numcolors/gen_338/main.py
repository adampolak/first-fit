# EVOLVE-BLOCK-START

def construct_intervals(seed_count=1):
  """
  Construct a sequence of intervals on the real line, in the order presented to FirstFit,
  to maximize FF colors divided by the clique number (omega).

  Interface preserved: construct_intervals(seed_count=1) -> list[(l, r)]
  """

  # Capacity guard to keep total intervals < 10000
  CAP = 9800

  # Four strong KT start templates (empirically best)
  TEMPLATE_BANK = [
    (2, 6, 10, 14),  # classic KT
    (1, 5, 9, 13),   # left-shifted
    (3, 7, 11, 15),  # right-shifted
    (4, 8, 12, 16),  # stretched-right
  ]

  # Seed with a single unit interval; multi-seed tends to inflate omega too early
  T = [(0, 1)]

  # ---------- Helpers ----------
  def _span_delta(current_T):
    lo = current_T[0][0]
    hi = current_T[0][1]
    # Single pass for speed
    for l, r in current_T:
      if l < lo: lo = l
      if r > hi: hi = r
    delta = hi - lo
    if delta <= 0:
      delta = 1
    return lo, hi, delta

  def _rebase_to_zero(current_T):
    if not current_T:
      return current_T
    lo, _, _ = _span_delta(current_T)
    if lo == 0:
      return current_T
    # Affine shift that preserves all intersections/order
    return [(l - lo, r - lo) for (l, r) in current_T]

  def _append_connectors(S, starts, delta):
    # Classic four connectors; preserves strong FF pressure with modest omega.
    s0, s1, s2, s3 = starts
    S.append(((s0 - 1) * delta, (s1 - 1) * delta))  # left cap
    S.append(((s2 + 2) * delta, (s3 + 2) * delta))  # right cap
    S.append(((s0 + 2) * delta, (s2 - 1) * delta))  # cross 1
    S.append(((s1 + 2) * delta, (s3 - 1) * delta))  # cross 2

  def _apply_round(current_T, starts, do_interleave=False, reverse_order=False):
    lo, hi, delta = _span_delta(current_T)

    # Build four translated blocks
    # Use list comprehension for speed
    base_vals = [s * delta - lo for s in starts]
    blocks = [[(l + base, r + base) for (l, r) in current_T] for base in base_vals]

    # Assemble with optional interleaving and reverse order to mix colors
    S = []
    if do_interleave:
      order = [0, 1, 2, 3]
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

    # Append classic connectors
    _append_connectors(S, starts, delta)
    return S

  def _cap_at(lo, span, a_frac, b_frac):
    L = lo + max(1, int(round(a_frac * span)))
    R = lo + max(1, int(round(b_frac * span)))
    if R <= L:
      R = L + 1
    return (L, R)

  # Predictive cap-aware number of KT rounds: n -> 4n + 4 per round
  def _max_rounds_within_cap(initial_size, max_rounds):
    sz = initial_size
    done = 0
    for _ in range(max(0, int(max_rounds))):
      nxt = 4 * sz + 4
      if nxt > CAP:
        break
      sz = nxt
      done += 1
    return done

  # ---------- Stage 1: KT spine with rotation and parity policies ----------
  target_rounds = _max_rounds_within_cap(len(T), 6)
  for ridx in range(target_rounds):
    # Rotate starts deterministically
    starts = TEMPLATE_BANK[ridx % len(TEMPLATE_BANK)]
    # Interleave on even rounds; reverse block order on odd rounds
    do_inter = (ridx % 2 == 0)
    rev = (ridx % 2 == 1)
    T = _apply_round(T, starts, do_interleave=do_inter, reverse_order=rev)
    if len(T) >= CAP:
      T = T[:CAP]
      return T
    # Keep numbers bounded without changing structure
    T = _rebase_to_zero(T)

  # Early exit if nearly at capacity
  if len(T) >= CAP - 8:
    return T

  # ---------- Stage 1.5: long-range cross-scale connectors ----------
  lo, hi, span = _span_delta(T)
  # Add cross-scale connectors: couple low and high regions to boost FF mixing
  cross_caps = [
    (lo + max(1, int(round(0.10 * span))), hi - max(1, int(round(0.10 * span)))),
    (lo + max(1, int(round(0.25 * span))), hi - max(1, int(round(0.25 * span)))),
    (lo + max(1, int(round(0.40 * span))), hi - max(1, int(round(0.40 * span)))),
    (lo + max(1, int(round(0.60 * span))), hi - max(1, int(round(0.60 * span)))),
  ]
  room = CAP - len(T)
  for (L, R) in cross_caps:
    if room <= 0:
      break
    if R <= L:
      R = L + 1
    T.append((L, R))
    room -= 1
  # Rebase to keep coordinates bounded
  T = _rebase_to_zero(T)

  # ---------- Micro-phase A: near-tail long caps ----------
  lo, hi, span = _span_delta(T)
  caps = [
    _cap_at(lo, span, 0.08, 0.60),
    _cap_at(lo, span, 0.25, 0.75),
    _cap_at(lo, span, 0.75, 0.92),
  ]

  # Insert near tail to couple many active colors without inflating omega
  def _insert_near_tail(seq, intervals):
    out = list(seq)
    for i, iv in enumerate(intervals):
      pos = len(out) - (i * 2 + 1)
      if pos < 0:
        out.append(iv)
      else:
        out.insert(pos, iv)
    return out

  room = CAP - len(T)
  if room > 0:
    T = _insert_near_tail(T, caps[:room])
    T = _rebase_to_zero(T)

  if len(T) >= CAP - 16:
    return T

  # ---------- Stage 2: fractional-span micro rounds (thin, even seed) ----------
  def _build_micro_delta_round(current_T, budget, iter_id=0, alt=False):
    if not current_T or budget <= 8:
      return []

    glo, ghi, G = _span_delta(current_T)

    # Thin seed from current_T, evenly spaced; slightly larger to improve mixing
    seed_sz = max(8, min(40, len(current_T) // 250))
    stride = max(1, len(current_T) // max(1, seed_sz))
    U = [current_T[i] for i in range(0, len(current_T), stride)][:seed_sz]
    if not U:
      return []

    ulo = min(l for l, _ in U)

    # Two families of windows: primary (A, with slight shifts) and alternate (B).
    if not alt:
      shift = (iter_id % 3) * 0.02
      window_fracs = [
        (0.12 + shift, 0.22 + shift),
        (0.35 + shift, 0.45 + shift),
        (0.58 + shift, 0.68 + shift),
        (0.80 + shift, 0.90 + shift),
      ]
      # Clamp windows conservatively inside the span
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
      tag = iter_id if not alt else (iter_id + 1)
      if ((int(round(fa * 100)) // 5) + tag) % 2 == 1:
        block = list(reversed(block))
      blocks.append(block)

    # Interleave micro-blocks (forward on even tag, reverse on odd tag)
    micro = []
    maxlen = max(len(b) for b in blocks)
    block_order = list(range(len(blocks)))
    tag = iter_id if not alt else (iter_id + 1)
    if tag % 2 == 1:
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
    if alt:
      # One longer cross reserved for alternate micro-phase
      micro_connectors.append((glo + int(round(0.18 * G)), glo + int(round(0.84 * G))))
    for a, b in micro_connectors:
      if b > a:
        micro.append((a, b))

    # Trim to available budget
    if len(micro) > budget:
      micro = micro[:budget]
    return micro

  # Execute up to two micro-rounds: primary then alternate
  room = CAP - len(T)
  if room > 8:
    microA = _build_micro_delta_round(T, room, iter_id=0, alt=False)
    if microA:
      if len(microA) > room:
        microA = microA[:room]
      T.extend(microA)
      T = _rebase_to_zero(T)

  room = CAP - len(T)
  if room > 8:
    microB = _build_micro_delta_round(T, room, iter_id=1, alt=True)
    if microB:
      if len(microB) > room:
        microB = microB[:room]
      T.extend(microB)
      T = _rebase_to_zero(T)

  # Micro-phase C: secondary fractional-window pass for finer coupling
  room = CAP - len(T)
  if room > 8:
    microC = _build_micro_delta_round(T, room, iter_id=2, alt=False)
    if microC:
      if len(microC) > room:
        microC = microC[:room]
      T.extend(microC)
      T = _rebase_to_zero(T)

  # Final capacity trim
  if len(T) > CAP:
    T = T[:CAP]

  return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()