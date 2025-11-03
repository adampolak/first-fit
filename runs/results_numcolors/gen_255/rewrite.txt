# EVOLVE-BLOCK-START

def construct_intervals(rounds=6,
                        rotate_starts=True,
                        reverse_block_parity=True,
                        interleave_blocks=True,
                        phase2_iters=2,
                        seed_count=1):
  """
  Deterministic KT-style spine + strengthened micro-phase construction.

  Signature kept identical to the previous program. Returns a list of integer
  open intervals (l, r), in the order they should be presented to FirstFit.
  """

  # Hard capacity guard (keep below 10000)
  CAP = 9800

  # Strong template bank (rotated deterministically)
  template_bank = [
    (2, 6, 10, 14),
    (1, 5, 9, 13),
    (3, 7, 11, 15),
    (4, 8, 12, 16),
  ]

  # Deterministic helper: small hash derived from integer inputs to avoid
  # use of global randomness while allowing reproducible per-round choices.
  def dseed(*args):
    h = 1469598103934665603
    for x in args:
      h ^= int(x) & 0xFFFFFFFFFFFFFFFF
      h = (h * 1099511628211) & 0xFFFFFFFFFFFFFFFF
    return h

  # Initialize seed set: single unit interval (classic KT) or a few sparse seeds
  if seed_count <= 1:
    T = [(0, 1)]
  else:
    seeds = min(4, max(1, int(seed_count)))
    step = 3
    T = [(i * step, i * step + 1) for i in range(seeds)]

  # Utility: span and integer delta (min delta = 1)
  def span_delta(seq):
    lo = min(l for l, r in seq)
    hi = max(r for l, r in seq)
    delta = hi - lo
    if delta <= 0:
      delta = 1
    return lo, hi, int(delta)

  # Classic KT connectors, shifted by lo to keep integer coordinates stable
  def append_kt_connectors(S, starts, delta, lo):
    s0, s1, s2, s3 = starts
    S.append(((s0 - 1) * delta + lo, (s1 - 1) * delta + lo))
    S.append(((s2 + 2) * delta + lo, (s3 + 2) * delta + lo))
    S.append(((s0 + 2) * delta + lo, (s2 - 1) * delta + lo))
    S.append(((s1 + 2) * delta + lo, (s3 - 1) * delta + lo))

  # Apply one KT-like round: translate blocks and append connectors.
  def apply_kt_round(current_T, starts, do_interleave=False, reverse_order=False, add_cross4=False):
    lo, hi, delta = span_delta(current_T)
    blocks = []
    for s in starts:
      base_off = s * delta - lo
      block = [(l + base_off, r + base_off) for (l, r) in current_T]
      blocks.append(block)

    S = []
    if do_interleave:
      order = list(range(4))
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

    # Classic connectors
    append_kt_connectors(S, starts, delta, lo)

    # Optional long cross4 connector to couple far layers (span 4..6)
    if add_cross4:
      s0 = starts[0]
      # produce a couple of long connectors spanning 4-6 delta windows
      A = (s0 + 4) * delta + lo
      B = (s0 + 6) * delta + lo
      if B > A:
        S.append((A, B))
      # symmetric one near the end
      s3 = starts[-1]
      A2 = (s3 - 6) * delta + lo
      B2 = (s3 - 4) * delta + lo
      if B2 > A2:
        S.append((A2, B2))

    return S

  # Predictive growth: size -> 4*size + 4 per KT round
  def max_rounds_within_cap(initial_size, max_rounds):
    sz = initial_size
    done = 0
    for _ in range(max(0, int(max_rounds))):
      nxt = 4 * sz + 4
      if nxt > CAP:
        break
      sz = nxt
      done += 1
    return done, sz

  # Stage 1: KT spine (deterministic rotation and interleaving policy)
  depth, _ = max_rounds_within_cap(len(T), rounds)
  for ridx in range(depth):
    starts = template_bank[ridx % len(template_bank)] if rotate_starts else template_bank[0]
    do_inter = bool(interleave_blocks and (ridx % 2 == 0))
    rev = bool(reverse_block_parity and (ridx % 2 == 1))
    # Occasionally add long-range connectors to improve cross-scale coupling (deterministic)
    add_cross4 = (ridx % 3 == 0)
    T = apply_kt_round(T, starts, do_interleave=do_inter, reverse_order=rev, add_cross4=add_cross4)
    if len(T) >= CAP:
      T = T[:CAP]
      return T

  if len(T) >= CAP - 20:
    return T[:CAP]

  # Insert a small deterministic set of early "blocker" pins to occupy low FF colors.
  # Each blocker is very short but placed at distinct fractional positions so they do not
  # create a large clique together (keeps omega small), yet they arrive early and
  # will tend to occupy small-numbered colors.
  def insert_blockers(seq, num_blockers=6):
    if num_blockers <= 0:
      return seq
    out = list(seq)
    lo, hi, delta = span_delta(seq)
    gap = max(1, delta // (num_blockers + 2))
    # place blockers near the left part to ensure they color early
    for i in range(num_blockers):
      a = lo + gap * (i + 1)
      # make blocker extremely short
      b = a + max(1, delta // 2000)
      # insert near the beginning (so they are presented early)
      out.insert(min(len(out), 1 + i * 2), (a, b))
    return out

  # Insert blockers only if we have room
  if len(T) < CAP - 40:
    T = insert_blockers(T, num_blockers=6)

  # Micro-phases: two different deterministic four-window micro-phases (CAP-gated).
  # Each micro-phase:
  #  - samples a thin seed U from current T
  #  - translates U into 4 windows; interleaves blocks
  #  - appends deterministic micro-connectors
  def thin_seed(current_T, max_seed):
    n = len(current_T)
    if n <= 0 or max_seed <= 0:
      return []
    step = max(1, n // max_seed)
    return current_T[::step][:max_seed]

  def build_micro_phase(current_T, budget, phase_id=0):
    if not current_T or budget <= 8:
      return []

    glo = min(l for l, r in current_T)
    ghi = max(r for l, r in current_T)
    G = max(1, ghi - glo)

    # deterministic seed size based on phase id and budget
    base_seed = max(8, min(36, (budget // 16) + (len(current_T) // 300)))
    seed_sz = base_seed if phase_id == 0 else max(8, min(40, (budget // 12) + (len(current_T) // 260)))
    U = thin_seed(current_T, seed_sz)
    if not U:
      return []

    ulo = min(l for l, r in U)

    # Two distinct window sets (phase 0 and phase 1), chosen deterministic and non-overlapping
    window_fracs = [
      [(0.12, 0.22), (0.32, 0.42), (0.58, 0.68), (0.78, 0.88)],  # phase 0
      [(0.05, 0.15), (0.28, 0.38), (0.60, 0.70), (0.82, 0.92)],  # phase 1 (shifted)
    ][phase_id % 2]

    # Slight deterministic micro-shift based on phase id to avoid perfect alignment
    micro_shift = ((dseed(len(current_T), phase_id) % 5) - 2) * 0.005

    blocks = []
    for (fa, fb) in window_fracs:
      fa_s = max(0.02, min(0.92, fa + micro_shift))
      win_lo = glo + int(round(fa_s * G))
      base = win_lo - ulo
      block = [(l + base, r + base) for (l, r) in U]
      # alternate reversal to break symmetry deterministically
      tag = int((fa * 100) // 5) + phase_id
      if tag % 2 == 1:
        block = list(reversed(block))
      blocks.append(block)

    # Interleave blocks to maximize color mixing
    micro = []
    maxlen = max(len(b) for b in blocks)
    block_order = list(range(len(blocks)))
    if (dseed(len(current_T), phase_id) & 1) == 1:
      block_order.reverse()
    for i in range(maxlen):
      for idx in block_order:
        blk = blocks[idx]
        if i < len(blk):
          micro.append(blk[i])

    # Deterministic connectors across windows (fractional spans)
    micro_connectors = [
      (glo + int(round(0.06 * G)), glo + int(round(0.28 * G))),
      (glo + int(round(0.34 * G)), glo + int(round(0.56 * G))),
      (glo + int(round(0.62 * G)), glo + int(round(0.92 * G))),
    ]
    for a, b in micro_connectors:
      if b > a:
        micro.append((a, b))

    # Add a long-range cross-scale connector (span ~4-6 delta) relative to the global span
    # This couples micro windows to the spine across scales but is placed only if it fits
    lo_sp, hi_sp, delta_sp = span_delta(current_T)
    # compute long connector endpoints as fractions of global span
    longA = lo_sp + max(1, int(round(0.14 * delta_sp)))
    longB = lo_sp + max(1, int(round(0.86 * delta_sp)))
    if longB > longA:
      micro.append((longA, longB))

    # trim to budget
    if len(micro) > budget:
      micro = micro[:budget]
    return micro

  # Execute up to phase2_iters micro rounds (but cap to 2)
  steps = min(max(0, int(phase2_iters)), 2)
  for pid in range(steps):
    room = CAP - len(T)
    if room <= 8:
      break
    micro = build_micro_phase(T, room, phase_id=pid)
    if not micro:
      break
    if len(micro) > room:
      micro = micro[:room]
    T.extend(micro)

  # Run a second distinct micro-phase set if there's room (this uses the other window set)
  room = CAP - len(T)
  if room > 8:
    micro2 = build_micro_phase(T, room, phase_id=steps)
    if micro2:
      if len(micro2) > room:
        micro2 = micro2[:room]
      T.extend(micro2)

  if len(T) >= CAP - 20:
    return T[:CAP]

  # Stage: add a few light-bridge translated blocks (smaller scale) to pressure FF further.
  lo, hi, delta = span_delta(T)
  if delta > 0:
    # pick a deterministic alternate start pattern
    starts_alt = (2, 4, 8, 12)
    # build one small translated pass with short bridges between consecutive blocks
    S = []
    for i, s in enumerate(starts_alt):
      base_off = s * delta - lo
      # Use a thin seed from current T to keep clique low
      U = thin_seed(T, max_seed=12)
      ulo = min(l for l, r in U) if U else lo
      block = [(l + (base_off - (ulo - lo)), r + (base_off - (ulo - lo))) for (l, r) in U] if U else []
      S.extend(block)
      # short bridge (tiny interval) connecting to next block to force cross-block overlaps
      if i + 1 < len(starts_alt):
        next_s = starts_alt[i + 1]
        a = int(lo + delta * (s + 0.5))
        b = int(lo + delta * (next_s + 0.5))
        # make bridge very short and integerized
        small_gap = max(1, delta // 1500)
        A = min(a, b) + small_gap
        B = max(a, b) - small_gap
        if B > A:
          S.append((A, B))

    # append a small handful of scaled connector gadgets
    for k in range(2):
      A = lo + int(round(0.18 * delta)) + k
      B = lo + int(round(0.82 * delta)) - k
      if B > A:
        S.append((A, B))

    # capacity guard
    room = CAP - len(T)
    if room > 0:
      if len(S) > room:
        S = S[:room]
      T.extend(S)

  # Final micro pins dispersed (very short intervals to tie windows)
  if len(T) < CAP - 4:
    glo2 = min(l for l, r in T)
    ghi2 = max(r for l, r in T)
    G2 = max(1, ghi2 - glo2)
    eps = max(1, G2 // 512)
    pins = []
    for frac in (0.09, 0.27, 0.46, 0.64, 0.83):
      mid = glo2 + int(round(frac * G2))
      pins.append((mid, mid + eps))
    room = CAP - len(T)
    if room > 0:
      pins = pins[:room]
      T.extend(pins)

  # Final trim just in case
  if len(T) > CAP:
    T = T[:CAP]
  return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()