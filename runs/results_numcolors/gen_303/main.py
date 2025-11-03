# EVOLVE-BLOCK-START

def construct_intervals():
  """
  Enhanced KT-style spine with deterministic interleaving, rotating templates,
  conservative CAP-aware densification, and two staged micro-phases.

  Returns:
    intervals: list of (l, r) integer tuples (open intervals) in presentation order.
  """

  # Capacity guard to keep total intervals < 10000
  CAP = 9800

  # Template bank (rotate among these to break symmetry)
  template_bank = [
    (2, 6, 10, 14),  # classic KT
    (1, 5, 9, 13),   # left shift
    (3, 7, 11, 15),  # right shift
    (4, 8, 12, 16),  # stretched-right
  ]

  # Deterministic base seed for per-round choices (used only for parity decisions)
  base_seed = 47

  # Seed with a single unit interval (keeps omega small early).
  T = [(0, 1)]

  # Helpers
  def _span(current):
    lo = min(l for l, r in current)
    hi = max(r for l, r in current)
    delta = hi - lo
    if delta <= 0:
      delta = 1
    return lo, hi, delta

  def thin_seed(current, max_seed):
    """Return a thin, evenly-spaced sample of current (size <= max_seed)."""
    n = len(current)
    if n == 0 or max_seed <= 0:
      return []
    step = max(1, n // max_seed)
    return current[::step][:max_seed]

  def apply_round(current_T, starts, ridx, densify=True):
    """Build the 4-block translated round with deterministic interleaving and connectors.
       Optionally add a small densify overlay to increase cross-block coupling."""
    lo, hi, delta = _span(current_T)

    # Build four translated blocks
    blocks = []
    for s in starts:
      base = s * delta - lo
      block = [(l + base, r + base) for (l, r) in current_T]
      blocks.append(block)

    # Deterministic interleaving: reverse order on odd parity to break symmetry.
    order = list(range(4))
    if ((base_seed + ridx) % 2) == 1:
      order.reverse()

    S = []
    maxlen = max(len(b) for b in blocks)
    for i in range(maxlen):
      for idx in order:
        blk = blocks[idx]
        if i < len(blk):
          S.append(blk[i])

    # Classic connectors (localized to the current delta scale)
    s0, s1, s2, s3 = starts
    connectors = [
      ((s0 - 1) * delta, (s1 - 1) * delta),  # left cap
      ((s2 + 2) * delta, (s3 + 2) * delta),  # right cap
      ((s0 + 2) * delta, (s2 - 1) * delta),  # cross 1
      ((s1 + 2) * delta, (s3 - 1) * delta),  # cross 2
    ]
    # Extra long-range connector family (conservative)
    long_cross = ((s0 - 2) * delta, (s3 + 3) * delta)
    connectors.append(long_cross)

    for a, b in connectors:
      if b > a:
        S.append((a, b))

    # Conservative densify overlay: small shifted copies of a thin seed.
    # This increases local coupling (FF pressure) while being CAP-aware and small.
    if densify and len(current_T) > 8:
      # Keep overlay size modest
      U = thin_seed(current_T, max_seed=min(24, max(8, len(current_T) // 300)))
      if U:
        shift = max(1, delta // 16)
        overlays = []
        # produce two small overlays with different shifts to spread interactions
        for mult in (1, 3):
          base_shift = mult * shift + (ridx % max(1, shift))
          for (l, r) in U:
            L = l + base_shift
            R = r + base_shift
            if R > L:
              overlays.append((L, R))
        # Insert overlays sparsely into S (interleaved positions)
        if overlays:
          step = max(3, len(S) // (len(overlays) + 1))
          pos = step
          for iv in overlays:
            if len(S) >= CAP:
              break
            if pos >= len(S):
              S.append(iv)
            else:
              S.insert(pos, iv)
            pos += step + 1

    return S

  # Stage 1: build KT spine, rotating templates, enforced interleaving every round
  max_rounds = 6
  for ridx in range(max_rounds):
    starts = template_bank[ridx % len(template_bank)]
    T = apply_round(T, starts, ridx, densify=True)
    if len(T) > CAP:
      break

  # If near capacity, return current sequence
  if len(T) >= CAP - 16:
    return T

  # Micro-phase A: a small tail of long-range caps to tie late colors while controlling omega
  lo, hi, delta = _span(T)
  def cap_at(a_frac, b_frac):
    L = lo + max(1, int(round(a_frac * delta)))
    R = lo + max(1, int(round(b_frac * delta)))
    if R <= L:
      R = L + 1
    return (L, R)
  caps = [cap_at(0.08, 0.60), cap_at(0.25, 0.75), cap_at(0.72, 0.95)]
  for c in caps:
    if len(T) >= CAP:
      break
    T.append(c)

  if len(T) >= CAP - 16:
    return T

  # Stage 2: two micro delta rounds (primary + alternate window families)
  def build_micro_round(current_T, budget, iter_id=0, alt=False):
    if not current_T or budget <= 8:
      return []

    glo = min(l for l, r in current_T)
    ghi = max(r for l, r in current_T)
    G = max(1, ghi - glo)

    # Choose seed size and windows depending on alt flag
    if not alt:
      seed_sz = max(8, min(32, len(current_T) // 300))
      shift = (iter_id % 3) * 0.02
      window_fracs = [
        (0.12 + shift, 0.22 + shift),
        (0.35 + shift, 0.45 + shift),
        (0.58 + shift, 0.68 + shift),
        (0.80 + shift, 0.90 + shift),
      ]
    else:
      seed_sz = max(8, min(40, len(current_T) // 250))
      window_fracs = [
        (0.05, 0.15),
        (0.28, 0.38),
        (0.60, 0.70),
        (0.82, 0.92),
      ]

    stride = max(1, len(current_T) // max(1, seed_sz))
    U = [current_T[i] for i in range(0, len(current_T), stride)][:seed_sz]
    if not U:
      return []

    ulo = min(l for l, r in U)
    # Slightly smaller delta2 for alt rounds to keep them focused
    delta2 = max(1, (G // 3) if not alt else (G // 2))

    # Build blocks and apply deterministic internal parity reversal
    blocks = []
    tag = iter_id + (1 if alt else 0)
    for (fa, fb) in window_fracs:
      win_lo = glo + int(round(fa * G))
      base = win_lo - ulo
      block = [(l + base, r + base) for (l, r) in U]
      if ((int(round(fa * 100)) // 5) + tag) % 2 == 1:
        block = list(reversed(block))
      blocks.append(block)

    # Interleave blocks (reverse order on odd tags)
    micro = []
    maxlen = max(len(b) for b in blocks)
    order = list(range(len(blocks)))
    if tag % 2 == 1:
      order.reverse()
    for i in range(maxlen):
      for idx in order:
        blk = blocks[idx]
        if i < len(blk):
          micro.append(blk[i])

    # Fractional-span connectors (micro-scale) and an optional longer cross for alternate
    micro_connectors = [
      (glo + int(round(0.08 * G)), glo + int(round(0.30 * G))),
      (glo + int(round(0.60 * G)), glo + int(round(0.92 * G))),
      (glo + int(round(0.26 * G)), glo + int(round(0.56 * G))),
      (glo + int(round(0.44 * G)), glo + int(round(0.78 * G))),
    ]
    if alt:
      micro_connectors.append((glo + int(round(0.18 * G)), glo + int(round(0.84 * G))))

    for a, b in micro_connectors:
      if b > a:
        micro.append((a, b))

    # Trim to budget
    if len(micro) > budget:
      micro = micro[:budget]
    return micro

  # Primary micro round
  room = CAP - len(T)
  if room > 8:
    micro1 = build_micro_round(T, max(0, room // 2), iter_id=0, alt=False)
    if micro1:
      avail = CAP - len(T)
      if len(micro1) > avail:
        micro1 = micro1[:avail]
      T.extend(micro1)

  # Alternate micro round
  room = CAP - len(T)
  if room > 8:
    micro2 = build_micro_round(T, max(0, room), iter_id=1, alt=True)
    if micro2:
      avail = CAP - len(T)
      if len(micro2) > avail:
        micro2 = micro2[:avail]
      T.extend(micro2)

  # Final trim to capacity
  if len(T) > CAP:
    T = T[:CAP]

  return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()