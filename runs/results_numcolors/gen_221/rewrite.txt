# EVOLVE-BLOCK-START

def construct_intervals(rounds=6,
                        rotate_starts=True,
                        reverse_block_parity=True,
                        interleave_blocks=True,
                        phase2_iters=2):
  """
  Duplex spine with cross4 connectors and two guarded micro-phases.
  Returns:
    intervals: list of (l, r) integer tuples, in FirstFit presentation order.
  """

  # Capacity guard: stay below 10000 with micro-phases
  CAP = 9800

  # Deterministic templates for the spine (rotated for robustness)
  template_bank = [
    (2, 6, 10, 14),  # classic KT
    (1, 5, 9, 13),   # left-shifted
    (3, 7, 11, 15),  # right-shifted
    (4, 8, 12, 16),  # stretched-right
  ]

  # Seed with a single unit interval to keep omega small initially
  T = [(0, 1)]

  # Helpers
  def _span(Tc):
    lo = min(l for l, r in Tc)
    hi = max(r for l, r in Tc)
    delta = hi - lo
    if delta <= 0:
      delta = 1
    return lo, hi, delta

  def _apply_spine_round(curr, starts, ridx, add_cross4=True):
    lo, hi, delta = _span(curr)
    # Duplex density: alternate K in {1, 2} to increase cross-scale coupling
    K = 1 if (ridx % 2 == 0) else 2

    # Build translated blocks
    blocks = []
    for s in starts:
      base = (s * K) * delta - lo
      blk = [(l + base, r + base) for (l, r) in curr]
      blocks.append(blk)

    # Parity-based assembly: even rounds interleave; odd rounds reversed sequential
    S = []
    if ridx % 2 == 0:
      maxlen = max(len(b) for b in blocks)
      for i in range(maxlen):
        for b in blocks:
          if i < len(b):
            S.append(b[i])
    else:
      for b in reversed(blocks):
        S.extend(b)

    # Classic connectors (conservative; scale respects K)
    s0, s1, s2, s3 = starts
    connectors = [
      ((s0 - 1) * K * delta, (s1 - 1) * K * delta),  # left cap
      ((s2 + 2) * K * delta, (s3 + 2) * K * delta),  # right cap
      ((s0 + 2) * K * delta, (s2 - 1) * K * delta),  # cross 1
      ((s1 + 2) * K * delta, (s3 - 1) * K * delta),  # cross 2
    ]
    for a, b in connectors:
      S.append((a, b))

    # Optional cross4 layer to couple distant spans; only on even rounds for safety
    if add_cross4 and (ridx % 2 == 0):
      # Two long bridges; lengths chosen to avoid creating large cliques
      c4 = [
        ((s0 + 4) * K * delta, (s3 + 4) * K * delta),
        ((s1 + 4) * K * delta, (s2 + 4) * K * delta),
      ]
      for a, b in c4:
        # Ensure valid orientation
        if b > a:
          S.append((a, b))
    return S

  # Plan how many rounds we can afford; estimate growth with connectors
  def est_next_size(sz, add_cross4):
    # Rough estimate: 4*sz + base connectors (4) + optional cross4 (2)
    return 4 * sz + (6 if add_cross4 else 4)

  # Stage 1: duplex spine with parity interleaving, capacity-guarded
  curr_size = len(T)
  for ridx in range(max(0, int(rounds))):
    starts = template_bank[ridx % len(template_bank)]
    add_c4 = True  # cross4 enabled but gated by CAP check below
    est = est_next_size(curr_size, add_c4)
    if est > CAP:
      # Try without cross4 if we are tight
      add_c4 = False
      est = est_next_size(curr_size, add_c4)
      if est > CAP:
        break
    T = _apply_spine_round(T, starts, ridx, add_cross4=add_c4)
    curr_size = len(T)
    if curr_size >= CAP - 32:
      break

  if len(T) >= CAP - 8:
    # Convert to integers and return early
    out = []
    for (l, r) in T[:CAP]:
      li = int(round(l))
      ri = int(round(r))
      if ri <= li:
        ri = li + 1
      out.append((li, ri))
    return out

  # Micro-phase seed helper: thin, evenly spaced, deterministic
  def thin_seed(seq, max_seed, div=300):
    n = len(seq)
    if n == 0 or max_seed <= 0:
      return []
    # Bound by both size-based and budget-based limits
    ss = max(8, min(max_seed, n // max(1, div)))
    step = max(1, n // ss)
    return seq[::step][:ss]

  # Micro-phase A: four-window translated blocks and connectors (delta2-like)
  def micro_phase_A(seq, budget, tag=0):
    if budget <= 8 or not seq:
      return []
    glo, ghi, G = _span(seq)
    U = thin_seed(seq, max_seed=40, div=250)
    if not U:
      return []
    ulo = min(l for l, r in U)
    # Window family A
    shift = (tag % 3) * 0.02
    windows = [
      (0.12 + shift, 0.22 + shift),
      (0.35 + shift, 0.45 + shift),
      (0.58 + shift, 0.68 + shift),
      (0.80 + shift, 0.90 + shift),
    ]
    windows = [(max(0.05, a), min(0.95, b)) for (a, b) in windows]

    blocks = []
    for idx, (fa, fb) in enumerate(windows):
      win_lo = glo + int(round(fa * G))
      base = win_lo - ulo
      block = [(l + base, r + base) for (l, r) in (U if idx % 2 == tag % 2 else list(reversed(U)))]
      blocks.append(block)

    # Interleave blocks with parity order to maximize FF coupling
    micro = []
    maxlen = max(len(b) for b in blocks)
    order = [0, 2, 1, 3] if (tag % 2 == 0) else [3, 1, 2, 0]
    for i in range(maxlen):
      for idx in order:
        blk = blocks[idx]
        if i < len(blk):
          micro.append(blk[i])

    # Fractional-span connectors (cross scale)
    connectors = [
      (glo + int(round(0.08 * G)), glo + int(round(0.30 * G))),
      (glo + int(round(0.60 * G)), glo + int(round(0.92 * G))),
      (glo + int(round(0.26 * G)), glo + int(round(0.56 * G))),
      (glo + int(round(0.44 * G)), glo + int(round(0.78 * G))),
      # cross4 long connector (guarded)
      (glo + int(round(0.18 * G)), glo + int(round(0.84 * G))),
    ]
    for a, b in connectors:
      if b > a:
        micro.append((a, b))

    return micro[:budget]

  # Micro-phase B: distinct windows with short "pin" intervals and light cross4 ties
  def micro_phase_B(seq, budget, tag=0):
    if budget <= 8 or not seq:
      return []
    glo, ghi, G = _span(seq)
    U = thin_seed(seq, max_seed=28, div=400)
    if not U:
      return []
    ulo = min(l for l, r in U)

    windows2 = [(0.05, 0.15), (0.28, 0.38), (0.60, 0.70), (0.82, 0.92)]
    eps = max(1, int(G // 512))  # very short pins

    blocks2 = []
    for widx, (fa, fb) in enumerate(windows2):
      win_lo = glo + int(round(fa * G))
      base = win_lo - ulo
      block = []
      # Staggered pins around midpoints of translated seed entries
      for idx, (l, r) in enumerate(U):
        mid = (l + r) // 2
        L = mid + base - (eps // 2) + (idx % 3)  # small stagger
        R = L + eps
        if R > L:
          block.append((L, R))
      # Optionally append a few full translated items (very sparse) to anchor coupling
      if widx % 2 == (tag % 2):
        for keep in range(0, len(U), max(4, len(U) // 8 or 1)):
          l, r = U[keep]
          block.append((l + base, r + base))
      blocks2.append(block)

    # Interleave with parity-based order
    micro2 = []
    order2 = [3, 1, 0, 2] if (tag % 2 == 0) else [2, 0, 1, 3]
    maxlen2 = max(len(b) for b in blocks2)
    for i in range(maxlen2):
      for idx in order2:
        blk = blocks2[idx]
        if i < len(blk):
          micro2.append(blk[i])

    # Short cross-window tie pins and a gentle long-range bridge
    ties = [
      (glo + int(round(0.09 * G)), glo + int(round(0.09 * G)) + eps),
      (glo + int(round(0.64 * G)), glo + int(round(0.64 * G)) + eps),
      (glo + int(round(0.22 * G)), glo + int(round(0.88 * G))),  # long gentle bridge
    ]
    for a, b in ties:
      if b > a:
        micro2.append((a, b))

    return micro2[:budget]

  # Apply micro phases; phase2_iters requests are respected up to 2
  room = CAP - len(T)
  if room > 8:
    mA = micro_phase_A(T, room, tag=0)
    if mA:
      if len(mA) > room:
        mA = mA[:room]
      T.extend(mA)

  room = CAP - len(T)
  if room > 8:
    # Second guarded micro-phase with distinct windows and pin shrinking
    mB = micro_phase_B(T, room, tag=1)
    if mB:
      if len(mB) > room:
        mB = mB[:room]
      T.extend(mB)

  # Final integer conversion (keep values modest and strictly increasing endpoints)
  out = []
  for (l, r) in T[:CAP]:
    li = int(round(l))
    ri = int(round(r))
    if ri <= li:
      ri = li + 1
    out.append((li, ri))
  return out

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()