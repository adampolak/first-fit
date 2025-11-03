# EVOLVE-BLOCK-START

def construct_intervals(rounds=6,
                        rotate_starts=True,
                        reverse_block_parity=True,
                        interleave_blocks=True,
                        phase2_iters=1):
  """
  New KT-style spine with six-template rotation, per-round deterministic densification,
  and a two-pass micro-phase. Keeps same signature/outputs as previous implementations.

  Parameters:
    Same external signature as original; internal behavior expanded as described in DESCRIPTION.
  Returns:
    List of (l, r) integer tuples (open intervals) in presentation order.
  """

  # Hard capacity guard
  CAP = 9800

  # Use a scaling factor so small fractional shifts aren't lost by integer truncation
  SCALE = 1000.0

  # Six deterministic start-templates to cycle through (more mixing than 4)
  template_bank = [
      (2, 6, 10, 14),
      (1, 5, 9, 13),
      (3, 7, 11, 15),
      (4, 8, 12, 16),
      (2, 5, 11, 14),   # slight variants to break symmetry
      (1, 6, 9, 16),
  ]

  # Seed with single unit interval (consistent with prior baselines)
  T = [(0.0, 1.0)]

  # Helpers to compute span and delta
  def _span_delta(current_T):
    lo = min(l for l, r in current_T)
    hi = max(r for l, r in current_T)
    delta = hi - lo
    if delta <= 0:
      delta = 1.0
    return lo, hi, delta

  # Append standard KT-style connectors for a given start-pattern and delta
  def _append_connectors(S, starts, delta):
    s0, s1, s2, s3 = starts
    # connectors chosen to link blocks but not to create huge pointwise overlap
    S.append(((s0 - 1) * delta, (s1 - 1) * delta))
    S.append(((s2 + 2) * delta, (s3 + 2) * delta))
    S.append(((s0 + 2) * delta, (s2 - 1) * delta))
    S.append(((s1 + 2) * delta, (s3 - 1) * delta))

  # Apply one spine round: create 4 translated blocks, assemble order, append connectors
  def apply_round(current_T, starts, do_interleave=False, reverse_order=False):
    lo, hi, delta = _span_delta(current_T)
    # Translate 4 blocks
    blocks = []
    for s in starts:
      base = s * delta - lo
      block = [(l + base, r + base) for (l, r) in current_T]
      blocks.append(block)

    # Assemble S either interleaving or sequential
    S = []
    if do_interleave:
      maxlen = max(len(b) for b in blocks)
      order = list(range(4))
      if reverse_order:
        order.reverse()
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

    # Append connectors (as short intervals on the same scale)
    _append_connectors(S, starts, delta)
    return S

  # Deterministic small densification step:
  # append a small fraction of intervals from the central block(s) shifted slightly
  # multiplier_frac controls fraction of block to duplicate (shifted); kept small to avoid omega spikes.
  def densify_after_round(S, fraction=0.08, jitter_scale=0.002):
    if not S:
      return S
    lo, hi, delta = _span_delta(S)
    # Estimate block size by chunk heuristic: try to find repeated connector patterns near end
    # If we can't identify blocks simply append nothing
    n = len(S)
    if n < 8:
      return S
    # Choose a central region inside (avoid the tail connectors)
    central_lo = lo + 0.25 * delta
    central_hi = lo + 0.75 * delta
    # Pick candidate intervals that lie mostly inside central region
    candidates = [iv for iv in S if (iv[0] >= central_lo and iv[1] <= central_hi)]
    if not candidates:
      return S
    add_count = max(1, int(len(candidates) * fraction))
    # Pick evenly spaced candidates to copy
    step = max(1, len(candidates) // add_count)
    picks = [candidates[i] for i in range(0, len(candidates), step)][:add_count]
    # Create shifted copies with tiny jitter; jitter relative to delta to avoid stacking at a point
    new_intervals = []
    base_jitter = jitter_scale * delta
    # deterministic sequence of small shifts to ensure reproducibility
    for idx, (l, r) in enumerate(picks):
      jitter = ((idx % 3) - 1) * base_jitter + (idx % 2) * (base_jitter / 3.0)
      nl = l + jitter
      nr = r + jitter
      if nr <= nl:
        nr = nl + 0.0001 * delta
      new_intervals.append((nl, nr))
    # Insert copies interleaved near the end (so they press on late FirstFit choices)
    out = list(S)
    insert_pos = max(0, len(out) - min(200, len(out)//10))
    for i, iv in enumerate(new_intervals):
      if len(out) >= CAP:
        break
      out.insert(insert_pos + i, iv)
    return out

  # Predictive rounds count (we will not exceed CAP)
  def round_next_size(sz):
    return 4 * sz + 4

  def max_rounds_within_cap(initial_size, max_rounds):
    sz = initial_size
    done = 0
    for _ in range(max(0, int(max_rounds))):
      nxt = round_next_size(sz)
      if nxt > CAP:
        break
      sz = nxt
      done += 1
    return done, sz

  # Stage 1: six-template spine (or up to rounds but capacity-guarded)
  depth, _ = max_rounds_within_cap(len(T), rounds)
  # Cap depth at 6 for the six-template rotation
  depth = min(depth, 6)
  for ridx in range(depth):
    # rotate among six templates if requested, else use the first template repeatedly
    starts = template_bank[ridx % len(template_bank)] if rotate_starts else template_bank[0]
    do_inter = bool(interleave_blocks and (ridx % 2 == 0))
    rev = bool(reverse_block_parity and (ridx % 2 == 1))
    T = apply_round(T, starts, do_interleave=do_inter, reverse_order=rev)

    # Per-round deterministic densification: small fraction to amplify FF pressure
    # We densify only if room remains comfortably under CAP
    if len(T) < CAP - 64:
      T = densify_after_round(T, fraction=0.08, jitter_scale=0.002)

    # Always apply a safety trim (keep floats for now; convert later)
    if len(T) >= CAP:
      break

  # If we are near capacity, early return
  if len(T) >= CAP - 8:
    # Normalize and integerize below
    pass
  else:
    # Micro-phase: two deterministic passes with different window families
    # Each pass uses a thin seed from T and interleaves translated blocks.
    def thin_seed(current_T, max_seed):
      if not current_T or max_seed <= 0:
        return []
      n = len(current_T)
      step = max(1, n // max_seed)
      return current_T[::step][:max_seed]

    def build_micro_pass(current_T, budget, pass_id=0):
      if not current_T or budget <= 8:
        return []
      glo = min(l for l, r in current_T)
      ghi = max(r for l, r in current_T)
      G = max(1.0, ghi - glo)

      # seed size slightly larger on second pass to improve coupling
      seed_sz = max(8, min(48, len(current_T) // (260 - 10 * pass_id)))
      U = thin_seed(current_T, seed_sz)
      if not U:
        return []

      ulo = min(l for l, r in U)

      # Two different window families (A/B) for pass_id 0 and (C/D) for pass_id 1
      if pass_id == 0:
        window_fracs = [
          (0.10, 0.20),
          (0.36, 0.46),
          (0.58, 0.68),
          (0.80, 0.90),
        ]
      else:
        window_fracs = [
          (0.05, 0.15),
          (0.28, 0.38),
          (0.50, 0.60),
          (0.70, 0.82),
        ]

      blocks = []
      for wi, (fa, fb) in enumerate(window_fracs):
        win_lo = glo + fa * G
        base = win_lo - ulo
        block = [(l + base, r + base) for (l, r) in U]
        # reverse every other block deterministically to break symmetry
        if wi % 2 == 1:
          block = list(reversed(block))
        blocks.append(block)

      # Interleave micro-blocks to maximize mixing
      micro = []
      maxlen = max(len(b) for b in blocks)
      order = list(range(len(blocks)))
      if pass_id % 2 == 1:
        order.reverse()
      for i in range(maxlen):
        for idx in order:
          blk = blocks[idx]
          if i < len(blk):
            micro.append(blk[i])

      # Fractional connectors at micro scale to tie windows
      connectors = [
        (glo + 0.07 * G, glo + 0.32 * G),
        (glo + 0.58 * G, glo + 0.92 * G),
        (glo + 0.24 * G, glo + 0.56 * G),
        (glo + 0.40 * G, glo + 0.78 * G),
      ]
      # Add one long-range cross connector on second pass
      if pass_id == 1:
        connectors.append((glo + 0.16 * G, glo + 0.84 * G))
      for a, b in connectors:
        if b > a:
          micro.append((a, b))

      # Cap micro length to budget
      if len(micro) > budget:
        micro = micro[:budget]
      return micro

    # Execute two micro passes, capacity-guarded
    micro_passes = 2
    for pid in range(micro_passes):
      room = CAP - len(T)
      if room <= 8:
        break
      micro = build_micro_pass(T, room, pass_id=pid)
      if not micro:
        break
      avail = CAP - len(T)
      if len(micro) > avail:
        micro = micro[:avail]
      # Insert micro pass blocks near the tail to pressure later FF choices
      insert_pos = max(0, len(T) - min(500, len(T)//6))
      for j, iv in enumerate(micro):
        if len(T) >= CAP:
          break
        T.insert(insert_pos + j, iv)

    # Append deterministic long-range connectors (4-6) computed from final span
    if len(T) < CAP - 6:
      lo, hi, delta = _span_delta(T)
      connectors_long = [
        (lo + 0.02 * delta, lo + 0.30 * delta),
        (lo + 0.24 * delta, lo + 0.52 * delta),
        (lo + 0.48 * delta, lo + 0.76 * delta),
        (lo + 0.70 * delta, lo + 0.98 * delta),
        # couple across thirds
        (lo + 0.12 * delta, lo + 0.88 * delta),
        (lo + 0.04 * delta, lo + 0.94 * delta),
      ]
      room = CAP - len(T)
      for iv in connectors_long[:room]:
        T.append(iv)

  # Final safety trim to capacity
  if len(T) > CAP:
    T = T[:CAP]

  # Normalize endpoints: shift to nonnegative and scale to integers preserving small shifts
  if not T:
    return []
  min_l = min(l for l, r in T)
  if min_l < 0:
    T = [(l - min_l, r - min_l) for (l, r) in T]

  # Scale and round into integer coordinates while ensuring r > l
  intervals = []
  for (l, r) in T:
    li = int(round(l * SCALE))
    ri = int(round(r * SCALE))
    if ri <= li:
      ri = li + 1
    # Keep coordinates reasonably bounded (they can be large; evaluator tolerated large values before)
    intervals.append((li, ri))
    if len(intervals) >= CAP:
      break

  return intervals

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()