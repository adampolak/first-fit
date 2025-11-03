# EVOLVE-BLOCK-START

def construct_intervals(seed_count=1, adaptive_micro=True):
  """
  six_template_backbone_two_pass_micro

  New deterministic generator:
    - six-template rotating backbone (rounds)
    - deterministic interleaving schedule (even rounds and multiples of 3)
    - per-round density multiplier to densify central blocks slightly
    - two-pass micro-phase (distinct window families)
    - deterministic long-range connectors appended after backbone

  Inputs maintained:
    seed_count (int): preserved parameter (keeps single-seed behaviour)
    adaptive_micro (bool): whether micro-phase uses multiple passes adaptively

  Output:
    list of (l, r) integer tuples in FirstFit presentation order
  """

  CAP = 9800  # hard cap on number of intervals

  # Six rotating templates (diverse start offsets)
  template_bank = [
    (2, 6, 10, 14),  # classic KT-like
    (1, 5, 9, 13),   # left-shifted
    (3, 7, 11, 15),  # right-shifted
    (4, 8, 12, 16),  # stretched-right
    (2, 4, 8, 12),   # compressed-left pair
    (3, 5, 9, 13),   # gentle left pack
  ]

  # Seed: single unit interval to permit many safe backbone rounds
  T = [(0, 1)] if seed_count <= 1 else [(0, 1)]

  # Helpers
  def _span(current_T):
    if not current_T:
      return 0, 0, 1
    lo = min(l for l, r in current_T)
    hi = max(r for l, r in current_T)
    delta = hi - lo
    if delta <= 0:
      delta = 1
    return lo, hi, delta

  def _even_sample(seq, k):
    """Evenly sample at most k items from seq (deterministic)."""
    if not seq or k <= 0:
      return []
    n = len(seq)
    step = max(1, n // k)
    return [seq[i] for i in range(0, n, step)][:k]

  def _append_classic_connectors(S, starts, delta):
    s0, s1, s2, s3 = starts
    # Four classic connectors (integer endpoints)
    S.append(((s0 - 1) * delta, (s1 - 1) * delta))
    S.append(((s2 + 2) * delta, (s3 + 2) * delta))
    S.append(((s0 + 2) * delta, (s2 - 1) * delta))
    S.append(((s1 + 2) * delta, (s3 - 1) * delta))

  def _apply_backbone_round(current_T, starts, do_interleave=False, reverse_order=False, densify=False, dens_mult=1.0):
    lo, hi, delta = _span(current_T)
    # Build 4 translated blocks from current_T using starts
    blocks = []
    for s in starts:
      base = s * delta - lo
      block = [(int(l + base), int(r + base)) for (l, r) in current_T]
      blocks.append(block)

    # Optional densification: duplicate a thin seed into central blocks with tiny shifts
    if densify and len(current_T) >= 8:
      # choose small seed and replicate into blocks 1 and 2
      seed_sz = min(24, max(8, len(current_T) // 200))
      seed = _even_sample(current_T, seed_sz)
      # a small integral shift to avoid exact duplicates
      shift_unit = max(1, delta // 1000)
      replicates = int(max(1, round((dens_mult - 1.0) * 10)))
      for idx in (1, 2):
        if idx < len(blocks):
          base = starts[idx] * delta - lo
          for rep in range(replicates):
            shift = rep * shift_unit
            extra = [(int(l + base + shift), int(r + base + shift)) for (l, r) in seed]
            # append extras interleaved at the end of that block to increase local pressure
            blocks[idx].extend(extra)

    # Assemble S with deterministic interleaving policy
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

    # Append classic connectors
    _append_classic_connectors(S, starts, delta)

    # Append one conservative extra cross connector tying distant blocks (keeps omega moderate)
    s0, s1, s2, s3 = starts
    # long cross intentionally kept outside immediate centers but still inside span scaling
    S.append((int((s0 + 3) * delta), int((s3 + 3) * delta)))

    return S

  # Predictive growth bound: size -> 4*size + 4 per backbone round
  def max_rounds_within_cap(initial_size, max_rounds):
    sz = initial_size
    done = 0
    for _ in range(max(0, int(max_rounds))):
      nxt = 4 * sz + 4
      if nxt > CAP:
        break
      sz = nxt
      done += 1
    return done

  # Choose backbone rounds: try up to 7 rounds but capacity-guarded
  MAX_BACKBONE_ROUNDS = 7
  depth = max_rounds_within_cap(len(T), MAX_BACKBONE_ROUNDS)
  for ridx in range(depth):
    # the six-template rotation
    starts = template_bank[ridx % len(template_bank)]
    # deterministic interleaving schedule: even rounds and multiples of 3 interleave
    do_interleave = (ridx % 2 == 0) or (ridx % 3 == 0)
    reverse_order = (ridx % 2 == 1)
    # per-round density multiplier (small, deterministic, bounded)
    dens_mult = 1.05 + 0.02 * ((ridx % 6) / 1.0)
    densify = (ridx >= 2)  # only densify after a couple of backbone rounds
    # capacity check
    if 4 * len(T) + 8 > CAP:
      break
    T = _apply_backbone_round(T, starts, do_interleave=do_interleave, reverse_order=reverse_order, densify=densify, dens_mult=dens_mult)

    # safety trim if grown beyond CAP
    if len(T) > CAP:
      T = T[:CAP]
      break

  # If nearly at capacity return
  if len(T) >= CAP - 16:
    return T[:CAP]

  # Deterministic long-range connectors appended after backbone (tie colors across scales)
  def append_long_connectors(current_T, max_add=6):
    if not current_T:
      return []
    lo, hi, delta = _span(current_T)
    span = max(1, hi - lo)
    # Spread connectors across the span: these are long intervals but chosen to avoid increasing clique maximum
    fracs = [
      (0.02, 0.36),
      (0.30, 0.62),
      (0.55, 0.88),
      (0.12, 0.78),
      (0.33, 0.66),
      (0.05, 0.95),
    ]
    connectors = []
    for (a, b) in fracs[:max_add]:
      L = lo + int(round(a * span))
      R = lo + int(round(b * span))
      if R <= L:
        R = L + 1
      connectors.append((L, R))
    return connectors

  # Append a small deterministic family of long connectors (capacity-aware)
  extra_connectors = append_long_connectors(T, max_add=4)
  for c in extra_connectors:
    if len(T) >= CAP:
      break
    T.append(c)

  if len(T) >= CAP - 16:
    return T[:CAP]

  # Two-pass micro-phase
  def build_micro_delta_round(current_T, budget, window_fracs, interleave=True, reverse_order=False, seed_limit=32, add_long_cross=False):
    if not current_T or budget <= 6:
      return []
    glo, ghi, G = _span(current_T)
    seed_sz = max(8, min(seed_limit, len(current_T) // 250))
    U = _even_sample(current_T, seed_sz)
    if not U:
      return []
    ulo = min(l for l, r in U)
    # Build micro-blocks aligned to windows
    blocks = []
    for (fa, fb) in window_fracs:
      win_lo = glo + int(round(max(0.0, fa) * G))
      base = win_lo - ulo
      block = [(int(l + base), int(r + base)) for (l, r) in U]
      blocks.append(block)
    # Interleave blocks deterministically
    micro = []
    if interleave:
      order = list(range(len(blocks)))
      if reverse_order:
        order.reverse()
      maxlen = max(len(b) for b in blocks)
      for i in range(maxlen):
        for idx in order:
          blk = blocks[idx]
          if i < len(blk):
            micro.append(blk[i])
    else:
      for blk in blocks:
        micro.extend(blk)
    # Add micro-scale connectors
    micro_connectors = [
      (glo + int(round(0.08 * G)), glo + int(round(0.30 * G))),
      (glo + int(round(0.60 * G)), glo + int(round(0.92 * G))),
      (glo + int(round(0.26 * G)), glo + int(round(0.56 * G))),
      (glo + int(round(0.44 * G)), glo + int(round(0.78 * G))),
    ]
    if add_long_cross:
      micro_connectors.append((glo + int(round(0.18 * G)), glo + int(round(0.84 * G))))
    for a, b in micro_connectors:
      if b > a:
        micro.append((a, b))
    if len(micro) > budget:
      micro = micro[:budget]
    return micro

  # Micro pass 1: primary windows
  lo, hi, delta = _span(T)
  Gspan = max(1, hi - lo)
  windows_A = [(0.12, 0.22), (0.35, 0.45), (0.58, 0.68), (0.80, 0.90)]
  room = CAP - len(T)
  if room > 8:
    microA = build_micro_delta_round(T, room, windows_A, interleave=True, reverse_order=False, seed_limit=40, add_long_cross=False)
    if microA:
      avail = CAP - len(T)
      if len(microA) > avail:
        microA = microA[:avail]
      T.extend(microA)

  # Micro pass 2: shifted windows, reverse interleaving + a long micro cross
  room = CAP - len(T)
  if room > 8:
    windows_B = [(0.05, 0.15), (0.28, 0.38), (0.60, 0.70), (0.82, 0.92)]
    microB = build_micro_delta_round(T, room, windows_B, interleave=True, reverse_order=True, seed_limit=48, add_long_cross=True)
    if microB:
      avail = CAP - len(T)
      if len(microB) > avail:
        microB = microB[:avail]
      T.extend(microB)

  # Adaptive micro extra pass if requested (small, capped)
  if adaptive_micro:
    room = CAP - len(T)
    if room > 8:
      # Another small micro pass with small shifts to catch more interactions
      windows_C = [(0.10, 0.18), (0.34, 0.42), (0.56, 0.66), (0.78, 0.88)]
      microC = build_micro_delta_round(T, room, windows_C, interleave=True, reverse_order=False, seed_limit=32, add_long_cross=False)
      if microC:
        avail = CAP - len(T)
        if len(microC) > avail:
          microC = microC[:avail]
        T.extend(microC)

  # Final append of a couple extra deterministic long connectors to tie remaining colors
  final_connectors = append_long_connectors(T, max_add=2)
  for c in final_connectors:
    if len(T) >= CAP:
      break
    T.append(c)

  # Final capacity trim
  if len(T) > CAP:
    T = T[:CAP]

  return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()