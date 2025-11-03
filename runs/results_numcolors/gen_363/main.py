# EVOLVE-BLOCK-START

def construct_intervals(seed_count=1):
  """
  Deterministic CAP-aware six-template backbone with two micro-passes and post-spine connectors.
  Returns a list of (l, r) integer open intervals, in FF arrival order.
  """
  CAP = 9800  # strict upper bound

  # Six strong start-template patterns to rotate across spine rounds.
  TEMPLATE_BANK = [
    (2, 6, 10, 14),  # classic KT
    (1, 5, 9, 13),   # left-shifted
    (3, 7, 11, 15),  # right-shifted
    (4, 8, 12, 16),  # stretched-right
    (2, 5, 11, 14),  # balanced inward
    (1, 7, 9, 15),   # mid-stretch
  ]

  # Deterministic base seed utilities (no randomness)
  BASE_SEED = 0xA5F1D3
  def lcg_next(x):
    return (x * 1103515245 + 12345) & 0x7FFFFFFF

  # Seed: maintain interface compatibility
  if seed_count <= 1:
    T = [(0, 1)]
  else:
    step = 3
    seeds = min(4, max(1, int(seed_count)))
    T = [(i * step, i * step + 1) for i in range(seeds)]

  # Utilities
  def span_delta(seq):
    lo = min(l for l, r in seq)
    hi = max(r for l, r in seq)
    d = hi - lo
    if d <= 0:
      d = 1
    return lo, hi, d

  def predict_next_size(sz):
    # KT-style recurrence (4 blocks + 4 connectors)
    return 4 * sz + 4

  def rounds_within_cap(initial_size, max_rounds):
    sz = initial_size
    done = 0
    for _ in range(int(max_rounds)):
      nxt = predict_next_size(sz)
      if nxt > CAP:
        break
      sz = nxt
      done += 1
    return done, sz

  def append_connectors(S, starts, delta, offset=0):
    s0, s1, s2, s3 = starts
    S.append((offset + (s0 - 1) * delta, offset + (s1 - 1) * delta))  # left cap
    S.append((offset + (s2 + 2) * delta, offset + (s3 + 2) * delta))  # right cap
    S.append((offset + (s0 + 2) * delta, offset + (s2 - 1) * delta))  # cross 1
    S.append((offset + (s1 + 2) * delta, offset + (s3 - 1) * delta))  # cross 2

  # Deterministic interleaving schedule: even rounds interleave, odd rounds sequential reversed.
  def apply_spine_round(current_T, starts, interleave, reverse_blocks, density_boost=False, boost_rate=0.08, seed_tag=0):
    lo, hi, delta = span_delta(current_T)
    # Build four translated blocks
    blocks = []
    for s in starts:
      base = s * delta - lo
      block = [(l + base, r + base) for (l, r) in current_T]
      blocks.append(block)

    # Build S
    S = []
    if interleave:
      maxlen = max(len(b) for b in blocks)
      order = list(range(len(blocks)))
      if reverse_blocks:
        order.reverse()
      for i in range(maxlen):
        for idx in order:
          blk = blocks[idx]
          if i < len(blk):
            S.append(blk[i])
    else:
      seq_blocks = list(reversed(blocks)) if reverse_blocks else blocks
      for blk in seq_blocks:
        S.extend(blk)

    # Per-round modest densification on central subset (CAP-aware)
    if density_boost and len(S) < CAP:
      # sample evenly spaced subset and re-embed them slightly left-shifted
      n = len(S)
      pick = max(8, min(64, n // 200))
      stride = max(1, n // max(1, pick))
      sample = [S[i] for i in range(n // 3, n - n // 3, stride)][:pick]
      shift = max(1, int(round(boost_rate * delta)))
      # Alternate direction deterministically by seed_tag
      shift = -shift if (seed_tag % 2 == 1) else shift
      densified = []
      for (l, r) in sample:
        L = l + shift
        R = r + shift
        if R > L:
          densified.append((L, R))
      # append densified block, trimming if near CAP
      free = CAP - len(S)
      if free > 0 and densified:
        if len(densified) > free:
          densified = densified[:free]
        S.extend(densified)

    # Append connectors at this scale
    append_connectors(S, starts, delta, offset=0)
    return S

  # Stage 1: Six-round backbone with six-template rotation and deterministic interleaving
  max_spine_rounds = 6
  depth, _ = rounds_within_cap(len(T), max_spine_rounds)
  seed = BASE_SEED
  for ridx in range(depth):
    starts = TEMPLATE_BANK[ridx % len(TEMPLATE_BANK)]
    interleave = (ridx % 2 == 0)
    reverse_blocks = (ridx % 2 == 1)
    density_boost = (ridx >= 1)  # enable from round 2 onwards
    T = apply_spine_round(T, starts, interleave, reverse_blocks, density_boost=density_boost, boost_rate=0.08, seed_tag=(seed ^ ridx))
    seed = lcg_next(seed)
    if len(T) >= CAP - 32:
      break

  # Stage 1.5: Post-spine deterministic long-range connectors (gate by CAP)
  if len(T) < CAP - 8:
    lo, hi, delta = span_delta(T)
    G = max(1, hi - lo)
    # 4–6 connectors spanning different fractions; avoid edges to bound omega
    LR = [
      (lo + int(0.07 * G), lo + int(0.29 * G)),
      (lo + int(0.21 * G), lo + int(0.51 * G)),
      (lo + int(0.49 * G), lo + int(0.79 * G)),
      (lo + int(0.68 * G), lo + int(0.91 * G)),
      (lo + int(0.13 * G), lo + int(0.87 * G)),  # one long cross
    ]
    # Insert near tail to maximize interaction with many active colors
    room = CAP - len(T)
    if room > 0:
      add = LR[:min(len(LR), room)]
      for i, iv in enumerate(add):
        pos = len(T) - (i * 3 + 1)
        if pos < 0:
          T.append(iv)
        else:
          T.insert(pos, iv)

  if len(T) >= CAP - 8:
    return _normalize_and_trim(T, CAP)

  # Helper: thin seed
  def thin_seed(current_T, max_seed):
    n = len(current_T)
    if n <= 0 or max_seed <= 0:
      return []
    step = max(1, n // max_seed)
    return current_T[::step][:max_seed]

  # Micro builder (window family, deterministic shifts)
  def build_micro(current_T, budget, iter_id, family_id):
    if budget <= 8 or not current_T:
      return []
    glo, ghi, G = span_delta(current_T)

    # seed selection tied to family_id
    seed_sz = max(12, min(40, len(current_T) // (250 if family_id == 1 else 300)))
    U = thin_seed(current_T, seed_sz)
    if not U:
      return []
    ulo = min(l for l, r in U)

    # Window families A and B with small deterministic shifts
    if family_id == 0:
      shift = ((iter_id + 1) % 3) * 0.015
      windows = [
        (0.12 + shift, 0.22 + shift),
        (0.34 + shift, 0.44 + shift),
        (0.56 + shift, 0.66 + shift),
        (0.78 + shift, 0.88 + shift),
      ]
    else:
      windows = [
        (0.08, 0.18),
        (0.28, 0.38),
        (0.60, 0.70),
        (0.82, 0.92),
      ]
    # clamp
    wf = []
    for a, b in windows:
      a = max(0.05, min(0.93, a))
      b = max(a + 0.01, min(0.95, b))
      wf.append((a, b))

    # Create translated micro-blocks
    blocks = []
    for (fa, fb) in wf:
      win_lo = glo + int(round(fa * G))
      base = win_lo - ulo
      block = [(l + base, r + base) for (l, r) in U]
      # symmetry-breaking reversal by parity
      if ((int(round(fa * 100)) // 5) + iter_id + family_id) % 2 == 1:
        block = list(reversed(block))
      blocks.append(block)

    # Interleave blocks in deterministic order
    micro = []
    maxlen = max(len(b) for b in blocks) if blocks else 0
    order = list(range(len(blocks)))
    if (iter_id + family_id) % 2 == 1:
      order.reverse()
    for i in range(maxlen):
      for idx in order:
        blk = blocks[idx]
        if i < len(blk):
          micro.append(blk[i])

    # Add fractional-span connectors for micro-level coupling
    connectors = [
      (glo + int(round(0.09 * G)), glo + int(round(0.30 * G))),
      (glo + int(round(0.28 * G)), glo + int(round(0.56 * G))),
      (glo + int(round(0.44 * G)), glo + int(round(0.78 * G))),
      (glo + int(round(0.62 * G)), glo + int(round(0.90 * G))),
    ]
    if family_id == 1:
      connectors.append((glo + int(round(0.18 * G)), glo + int(round(0.84 * G))))
    for a, b in connectors:
      if b > a:
        micro.append((a, b))

    # Trim to budget
    if len(micro) > budget:
      micro = micro[:budget]
    return micro

  # Stage 2: Two deterministic micro passes with separate window families
  remaining = CAP - len(T)
  if remaining > 12:
    # Pass A
    budgetA = max(8, remaining // 3)
    microA = build_micro(T, budgetA, iter_id=0, family_id=0)
    if microA:
      free = CAP - len(T)
      if len(microA) > free:
        microA = microA[:free]
      T.extend(microA)

  remaining = CAP - len(T)
  if remaining > 12:
    # Pass B
    budgetB = max(8, remaining // 2)
    microB = build_micro(T, budgetB, iter_id=1, family_id=1)
    if microB:
      free = CAP - len(T)
      if len(microB) > free:
        microB = microB[:free]
      T.extend(microB)

  # Final deterministic interleaved tail caps to tie distant regions lightly
  if len(T) < CAP:
    lo, hi, d = span_delta(T)
    d2 = max(1, d // 4)
    tail = [
      (lo + d2, lo + 3 * d2),
      (lo + int(0.35 * d), lo + int(0.72 * d)),
      (hi - 3 * d2, hi - d2),
    ]
    room = CAP - len(T)
    add = tail[:room]
    for i, iv in enumerate(add):
      pos = len(T) - (2 * i + 1)
      if pos < 0:
        T.append(iv)
      else:
        T.insert(pos, iv)

  return _normalize_and_trim(T, CAP)


def _normalize_and_trim(T, CAP):
  # Normalize to non-negative integers and ensure r > l
  if not T:
    return []
  min_l = min(l for l, r in T)
  if min_l < 0:
    T = [(l - min_l, r - min_l) for (l, r) in T]
  out = []
  for (l, r) in T:
    li = int(l)
    ri = int(r)
    if ri <= li:
      ri = li + 1
    out.append((li, ri))
  if len(out) > CAP:
    out = out[:CAP]
  return out

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()