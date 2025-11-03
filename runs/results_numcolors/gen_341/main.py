# EVOLVE-BLOCK-START

def construct_intervals(seed_count=1):
  """
  Construct a sequence of intervals of the real line, in FirstFit presentation order,
  to maximize FF colors divided by the clique number (omega). Returns:
    list[(l, r)] where each (l, r) is an open interval with integer endpoints.
  """

  CAP = 9800  # hard capacity guard to keep total intervals under 10000

  # Six-template backbone (diversified KT starts)
  TEMPLATES_6 = [
    (2, 6, 10, 14),  # classic KT
    (1, 5, 9, 13),   # left-shifted
    (3, 7, 11, 15),  # right-shifted
    (4, 8, 12, 16),  # stretched-right
    (2, 4, 8, 12),   # compressed left pair
    (2, 8, 10, 12),  # inner symmetric
  ]

  # Seed: single unit interval to avoid early omega inflation.
  # Multi-seed is supported but defaults to one seed for best behavior.
  if seed_count <= 1:
    T = [(0, 1)]
  else:
    step = 3
    T = [(i * step, i * step + 1) for i in range(max(1, int(seed_count)))]

  # Helpers
  def _span(seq):
    lo = min(l for l, _ in seq)
    hi = max(r for _, r in seq)
    d = hi - lo
    return lo, hi, (d if d > 0 else 1)

  def _append_connectors(S, starts, delta):
    # Classic four connectors shown to keep strong FF pressure and modest omega.
    s0, s1, s2, s3 = starts
    S.append(((s0 - 1) * delta, (s1 - 1) * delta))  # left cap
    S.append(((s2 + 2) * delta, (s3 + 2) * delta))  # right cap
    S.append(((s0 + 2) * delta, (s2 - 1) * delta))  # cross 1
    S.append(((s1 + 2) * delta, (s3 - 1) * delta))  # cross 2)

  def _apply_round(cur, starts, do_interleave=False, reverse_order=False):
    lo, hi, delta = _span(cur)
    # Build four translated blocks
    blocks = []
    for s in starts:
      base = s * delta - lo
      blk = [(l + base, r + base) for (l, r) in cur]
      blocks.append(blk)
    # Assemble presentation order
    S = []
    if do_interleave:
      maxlen = max(len(b) for b in blocks)
      order = list(range(4))
      if reverse_order:
        order.reverse()
      for i in range(maxlen):
        for idx in order:
          b = blocks[idx]
          if i < len(b):
            S.append(b[i])
    else:
      if reverse_order:
        blocks = list(reversed(blocks))
      for b in blocks:
        S.extend(b)
    _append_connectors(S, starts, delta)
    return S

  def _template_for_round(ridx):
    # Deterministic choice among six templates (no randomness, stable across runs)
    # Lightweight round-based indexing for template diversity
    idx = (1729 * (ridx + 1) + 2654435761) % 6
    return TEMPLATES_6[idx]

  # Stage 1: Six KT-style rounds with controlled interleaving
  # Size after k rounds with seed 1: n_k = 4^k + 4(4^{k-1}+...+1) = 4n + 4 recurrence; after k=6 it's 9556.
  for ridx in range(6):
    # Ensure next round fits comfortably (exact formula: 4*len(T) + 4)
    if 4 * len(T) + 4 > CAP:
      break
    starts = _template_for_round(ridx)
    do_interleave = (ridx % 2 == 0)   # interleave on even rounds
    reverse_order = (ridx % 2 == 1)   # reverse block order on odd rounds
    T = _apply_round(T, starts, do_interleave=do_interleave, reverse_order=reverse_order)
    if len(T) >= CAP:
      return T[:CAP]

  # Early exit if near capacity
  if len(T) >= CAP - 8:
    return T[:CAP]

  # Deterministic long-range caps after completing the spine (pairwise disjoint)
  # These are placed to increase coupling without inflating omega significantly.
  lo, hi, delta = _span(T)
  span = max(1, hi - lo)
  def frac_iv(a, b):
    L = lo + max(1, int(round(a * span)))
    R = lo + max(1, int(round(b * span)))
    return (L, R) if R > L else (L, L + 1)

  spine_caps = [
    frac_iv(0.06, 0.18),
    frac_iv(0.28, 0.40),
    frac_iv(0.46, 0.58),
    frac_iv(0.64, 0.76),
  ]
  for c in spine_caps:
    if len(T) >= CAP:
      break
    T.append(c)

  if len(T) >= CAP - 16:
    return T[:CAP]

  # Thin, evenly spaced seed from current T
  def thin_seed(cur, budget_cap):
    n = len(cur)
    if n == 0 or budget_cap <= 0:
      return []
    want = max(8, min(budget_cap, 40))
    step = max(1, n // want)
    return cur[::step][:want]

  # Build a micro-phase wave anchored to fractional windows across the span.
  def build_micro_wave(current_T, budget, windows, add_connectors=True, alt_tag=0):
    if not current_T or budget <= 8:
      return []

    glo = min(l for l, _ in current_T)
    ghi = max(r for _, r in current_T)
    G = max(1, ghi - glo)

    # Thin, evenly spaced seed (bounded to keep omega modest)
    seed = thin_seed(current_T, min(40, max(8, len(current_T) // 250)))
    if not seed:
      return []
    ulo = min(l for l, _ in seed)

    # Build translated micro-blocks at given windows; use staggered reversal to break symmetry
    blocks = []
    for j, (fa, fb) in enumerate(windows):
      win_lo = glo + int(round(fa * G))
      base = win_lo - ulo
      blk = [(l + base, r + base) for (l, r) in seed]
      if ((j + alt_tag) % 2) == 1:
        blk.reverse()
      blocks.append(blk)

    # Interleave blocks; alternate order with alt_tag for diversity
    micro = []
    maxlen = max(len(b) for b in blocks)
    order = list(range(len(blocks)))
    if alt_tag % 2 == 1:
      order.reverse()
    for i in range(maxlen):
      for idx in order:
        b = blocks[idx]
        if i < len(b):
          micro.append(b[i])

    if add_connectors:
      # Fractional-span connectors; moderate overlap to raise FF but control omega.
      cons = [
        (glo + int(round(0.08 * G)), glo + int(round(0.30 * G))),  # left cap
        (glo + int(round(0.26 * G)), glo + int(round(0.56 * G))),  # mid cross
        (glo + int(round(0.44 * G)), glo + int(round(0.78 * G))),  # long cross
        (glo + int(round(0.60 * G)), glo + int(round(0.92 * G))),  # right cap
      ]
      for a, b in cons:
        if b > a:
          micro.append((a, b))

    return micro[:budget]

  # Two-phase micro budget: split remaining capacity deterministically
  remaining = CAP - len(T)
  if remaining <= 0:
    return T[:CAP]

  budget1 = max(0, int(0.6 * remaining))
  budget2 = max(0, remaining - budget1)

  # Phase A windows (inside windows; avoid boundary to control omega)
  windows_A = [(0.12, 0.22), (0.35, 0.45), (0.58, 0.68), (0.80, 0.90)]
  # Phase B windows (shifted pattern)
  windows_B = [(0.05, 0.15), (0.28, 0.38), (0.60, 0.70), (0.82, 0.92)]

  waveA = build_micro_wave(T, budget1, windows_A, add_connectors=True, alt_tag=0) if budget1 > 0 else []
  if waveA:
    T.extend(waveA[:max(0, CAP - len(T))])

  remaining = CAP - len(T)
  if remaining > 0 and budget2 > 0:
    waveB = build_micro_wave(T, remaining, windows_B, add_connectors=True, alt_tag=1)
    if waveB:
      T.extend(waveB[:max(0, CAP - len(T))])

  # Final capacity clamp
  if len(T) > CAP:
    T = T[:CAP]
  return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()