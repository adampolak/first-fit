# EVOLVE-BLOCK-START

def _span_delta(T):
  lo = min(l for l, _ in T)
  hi = max(r for _, r in T)
  delta = hi - lo
  if delta <= 0:
    delta = 1
  return lo, hi, delta

def _build_blocks_sequential(T, starts, delta, lo):
  """
  Build four translated blocks sequentially (no interleaving).
  This ordering is known to amplify FirstFit pressure while keeping omega modest.
  """
  S = []
  for s in starts:
    base = s * delta - lo
    S.extend((l + base, r + base) for (l, r) in T)
  return S

def _append_connectors(S, starts, delta):
  """
  Append the classic four connectors used in KT-style constructions.
  """
  s0, s1, s2, s3 = starts
  S.append(((s0 - 1) * delta, (s1 - 1) * delta))  # left cap
  S.append(((s2 + 2) * delta, (s3 + 2) * delta))  # right cap
  S.append(((s0 + 2) * delta, (s2 - 1) * delta))  # cross 1
  S.append(((s1 + 2) * delta, (s3 - 1) * delta))  # cross 2

def _thin_seed_with_offset(T, target_sz, offset):
  """
  Evenly spaced thin seed selection of size up to target_sz, with a deterministic
  offset to diversify the chosen subset across micro-phases.
  """
  n = len(T)
  if n == 0 or target_sz <= 0:
    return []
  step = max(1, n // target_sz)
  # Clamp offset into [0, step-1]
  if step > 0:
    off = offset % step
  else:
    off = 0
  U = []
  idx = off
  while len(U) < target_sz and idx < n:
    U.append(T[idx])
    idx += step
  return U

def _micro_round(current_T, windows, parity_reverse=False, seed_tag=0, budget=0):
  """
  Build a micro-phase by taking a thin seed U and translating it into the given windows.
  Interleave blocks to increase FF pressure, add fractional-span connectors including cross4.

  Parameters:
    current_T: existing interval sequence
    windows: list of (a_frac, b_frac) pairs, each 0 < a < b < 1
    parity_reverse: if True, reverse the block interleaving order
    seed_tag: small integer to derive a deterministic offset for seed selection
    budget: maximum number of intervals to emit (0 => unbounded; will be capped by caller)
  """
  if not current_T:
    return []

  glo, ghi, G = _span_delta(current_T)

  # Thin seed size adapted to structure size; keep small to protect omega
  target_sz = max(8, min(32, len(current_T) // 300))
  # Deterministic derived seed offset (avoid RNG; use multiplicative hash)
  BASE = 1337
  M = 2654435761  # Knuth's multiplicative constant
  hashed = (BASE * M + seed_tag) & 0xFFFFFFFF
  U = _thin_seed_with_offset(current_T, target_sz, offset=hashed)
  if not U:
    return []

  ulo = min(l for l, _ in U)

  # Build translated blocks, aligned to window left endpoints
  blocks = []
  for (fa, fb) in windows:
    fa = max(0.05, min(0.90, fa))
    fb = max(0.10, min(0.95, fb))
    if fb <= fa:
      fb = fa + 0.05
    win_lo = glo + int(round(fa * G))
    base = win_lo - ulo
    block = [(l + base, r + base) for (l, r) in U]
    blocks.append(block)

  # Interleave micro-blocks with optional reversed block order
  interleaved = []
  order = list(range(len(blocks)))
  if parity_reverse:
    order.reverse()
  maxlen = max(len(b) for b in blocks)
  for i in range(maxlen):
    for idx in order:
      b = blocks[idx]
      if i < len(b):
        interleaved.append(b[i])

  # Deterministic fractional-span connectors, including two long-range cross4 connectors
  con = [
    (glo + int(round(0.08 * G)), glo + int(round(0.30 * G))),  # left cap
    (glo + int(round(0.60 * G)), glo + int(round(0.92 * G))),  # right cap
    (glo + int(round(0.26 * G)), glo + int(round(0.56 * G))),  # cross 1
    (glo + int(round(0.44 * G)), glo + int(round(0.78 * G))),  # cross 2
    (glo + int(round(0.18 * G)), glo + int(round(0.84 * G))),  # cross4a
    (glo + int(round(0.12 * G)), glo + int(round(0.88 * G))),  # cross4b
  ]
  for a, b in con:
    if b > a:
      interleaved.append((a, b))

  if budget and len(interleaved) > budget:
    interleaved = interleaved[:budget]
  return interleaved

def construct_intervals(seed_count=1):
  """
  Construct a sequence of open intervals (l, r) in FirstFit arrival order,
  aiming to maximize FirstFit colors relative to the clique number.

  Returns: list of (l, r) tuples with integer endpoints.
  """
  CAP = 9800

  # Stage 1: six-round KT-style spine with rotating templates, sequential blocks, classic connectors.
  template_bank = [
    (2, 6, 10, 14),  # classic KT
    (1, 5, 9, 13),   # left-shifted
    (3, 7, 11, 15),  # right-shifted
    (4, 8, 12, 16),  # stretched right
  ]

  # Seed: a single unit interval to avoid early omega inflation.
  T = [(0, 1)]

  # Build the spine
  rounds = 6
  for ridx in range(rounds):
    lo, hi, delta = _span_delta(T)
    starts = template_bank[ridx % len(template_bank)]
    S = _build_blocks_sequential(T, starts, delta, lo)
    _append_connectors(S, starts, delta)
    T = S
    if len(T) >= CAP:
      # Safety guard (should not trigger for the classic recurrence)
      T = T[:CAP]
      return T

  # Stage 2: dual micro-phase with distinct window sets and parity interleaving.
  if len(T) < CAP - 8:
    glo, ghi, G = _span_delta(T)
    room = CAP - len(T)

    # Micro A: interior windows; forward interleaving
    windows_A = [(0.12, 0.22), (0.35, 0.45), (0.58, 0.68), (0.80, 0.90)]
    microA_budget = max(0, room // 2)
    microA = _micro_round(T, windows_A, parity_reverse=False, seed_tag=1, budget=microA_budget)
    if microA:
      addA = min(len(microA), CAP - len(T))
      T.extend(microA[:addA])

  # Micro B: distinct windows; reversed interleaving; includes cross4 connectors
  if len(T) < CAP - 4:
    room = CAP - len(T)
    windows_B = [(0.05, 0.15), (0.28, 0.38), (0.60, 0.70), (0.82, 0.92)]
    microB = _micro_round(T, windows_B, parity_reverse=True, seed_tag=2, budget=room)
    if microB:
      addB = min(len(microB), CAP - len(T))
      T.extend(microB[:addB])

  # Hard cap
  if len(T) > CAP:
    T = T[:CAP]
  return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()