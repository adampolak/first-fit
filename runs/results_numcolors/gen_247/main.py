# EVOLVE-BLOCK-START

def construct_intervals(seed_count=1):
  """
  Construct a sequence of open intervals (l, r) for FirstFit, aiming to maximize
  FF colors divided by the clique number (<= 10 target). Deterministic and CAP-gated.
  """

  # Hard capacity guard
  CAP = 9800
  # Margin that must remain free for any micro-phase to run
  CAP_MARGIN = 32

  # Global deterministic seed and a tiny LCG for reproducible shifts
  BASE_SEED = 91138233

  class DRand:
    def __init__(self, s):
      self.state = (s ^ 0x9E3779B97F4A7C15) & ((1 << 64) - 1)
    def next(self):
      # LCG parameters (PCG-style constants)
      self.state = (self.state * 6364136223846793005 + 1442695040888963407) & ((1 << 64) - 1)
      return self.state
    def rand01(self):
      return ((self.next() >> 11) & ((1 << 53) - 1)) / float(1 << 53)

  # Rotating four-start templates for the KT-like spine
  template_bank = [
    (2, 6, 10, 14),  # classic KT
    (1, 5, 9, 13),   # left-shifted
    (3, 7, 11, 15),  # right-shifted
    (4, 8, 12, 16),  # stretched-right
  ]

  # Seed: single short interval to keep initial omega minimal
  T = [(0, 1)]

  # Helpers
  def _span(current_T):
    lo = min(l for l, r in current_T)
    hi = max(r for l, r in current_T)
    delta = hi - lo
    if delta <= 0:
      delta = 1
    return lo, hi, delta

  def _apply_round(current_T, starts, do_interleave=False, reverse_order=False):
    lo, hi, delta = _span(current_T)

    # Build four translated blocks
    blocks = []
    for s in starts:
      base = s * delta - lo
      block = [(l + base, r + base) for (l, r) in current_T]
      blocks.append(block)

    # Assemble with interleaving policy
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

    # Classic four connectors (safe, proven pattern)
    s0, s1, s2, s3 = starts
    S.append(((s0 - 1) * delta, (s1 - 1) * delta))  # left cap
    S.append(((s2 + 2) * delta, (s3 + 2) * delta))  # right cap
    S.append(((s0 + 2) * delta, (s2 - 1) * delta))  # cross 1
    S.append(((s1 + 2) * delta, (s3 - 1) * delta))  # cross 2
    return S

  # Stage 1: CAP-disciplined spine with parity-based interleaving
  # Predictive size: sz -> 4*sz + 4 per round
  def _can_apply_round(sz):
    return (4 * sz + 4) <= CAP

  ridx = 0
  rng = DRand(BASE_SEED ^ 0xA5A5A5A5)
  while ridx < 6 and _can_apply_round(len(T)):
    starts = template_bank[ridx % len(template_bank)]
    do_inter = (ridx % 2 == 0)
    rev = (ridx % 2 == 1)
    T = _apply_round(T, starts, do_interleave=do_inter, reverse_order=rev)
    ridx += 1
    if len(T) >= CAP:
      T = T[:CAP]
      return T

  # Early exit if nearly at capacity
  if len(T) >= CAP - CAP_MARGIN:
    return T

  # Connector layer A: deterministic long-range fractional caps across the global span
  # These "caps" and "crosses" are computed from the entire span to couple many active colors.
  def add_fractional_connectors(seq, max_add=6):
    lo, hi, delta = _span(seq)
    G = hi - lo
    if G <= 0:
      return seq
    # Deterministic small shifts from BASE_SEED to avoid exact alignment with spine joints
    r = DRand(BASE_SEED ^ 0xBEEF)
    fracs = [
      (0.08 + 0.01 * (r.rand01() * 0.5), 0.30 + 0.01 * (r.rand01() * 0.5)),
      (0.26 + 0.01 * (r.rand01() * 0.5), 0.56 + 0.01 * (r.rand01() * 0.5)),
      (0.44 + 0.01 * (r.rand01() * 0.5), 0.78 + 0.01 * (r.rand01() * 0.5)),
      (0.62 + 0.01 * (r.rand01() * 0.5), 0.90 + 0.01 * (r.rand01() * 0.5)),
      (0.18 + 0.01 * (r.rand01() * 0.5), 0.84 + 0.01 * (r.rand01() * 0.5)),  # long cross
      (0.33 + 0.01 * (r.rand01() * 0.5), 0.67 + 0.01 * (r.rand01() * 0.5)),  # mid cross
    ]
    added = 0
    for (a, b) in fracs:
      L = lo + int(round(max(0.02, min(0.93, a)) * G))
      R = lo + int(round(max(0.07, min(0.98, b)) * G))
      if R <= L:
        R = L + 1
      if len(seq) < CAP and added < max_add:
        seq.append((L, R))
        added += 1
      else:
        break
    return seq

  if len(T) <= CAP - CAP_MARGIN:
    T = add_fractional_connectors(T, max_add=4)

  if len(T) >= CAP - CAP_MARGIN:
    return T

  # Micro-phase: thin, deterministic four-window micro-round
  def build_micro_round(current_T, budget, seed_shift=0, windows=None):
    if not current_T or budget <= 8:
      return []

    glo = min(l for l, r in current_T)
    ghi = max(r for l, r in current_T)
    G = max(1, ghi - glo)

    # Thin, evenly spaced seed (slightly larger bound than prior code)
    max_seed = max(12, min(48, len(current_T) // 240))
    stride = max(1, len(current_T) // max(1, max_seed))
    U = [current_T[i] for i in range(0, len(current_T), stride)][:max_seed]
    if not U:
      return []

    ulo = min(l for l, r in U)

    # Deterministic window shifts from the global seed chain
    r = DRand(BASE_SEED ^ (0x1234ABCD + seed_shift))
    # Default four-window family if none provided
    if windows is None:
      windows = [(0.12, 0.22), (0.35, 0.45), (0.58, 0.68), (0.80, 0.90)]
    # Apply tiny deterministic jitter to avoid resonance
    jittered = []
    for (fa, fb) in windows:
      j = (r.rand01() - 0.5) * 0.02  # +/- 0.01 jitter
      a = max(0.05, min(0.92, fa + j))
      b = max(a + 0.05, min(0.97, fb + j))
      jittered.append((a, b))

    # Build translated micro-blocks aligned to these windows
    blocks = []
    for (fa, fb) in jittered:
      win_lo = glo + int(round(fa * G))
      base = win_lo - ulo
      block = [(l + base, r + base) for (l, r) in U]
      # Reverse every other block deterministically to diversify FF’s local history
      if int(round(fa * 100)) % 2 == 1:
        block = list(reversed(block))
      blocks.append(block)

    # Interleave with round-dependent permutation
    micro = []
    order = list(range(len(blocks)))
    if seed_shift % 2 == 1:
      order.reverse()
    maxlen = max(len(b) for b in blocks)
    for i in range(maxlen):
      for idx in order:
        blk = blocks[idx]
        if i < len(blk):
          micro.append(blk[i])

    # Fractional connectors at micro scale (conservative set)
    connectors = [
      (glo + int(round(0.08 * G)), glo + int(round(0.30 * G))),
      (glo + int(round(0.26 * G)), glo + int(round(0.56 * G))),
      (glo + int(round(0.44 * G)), glo + int(round(0.78 * G))),
      (glo + int(round(0.60 * G)), glo + int(round(0.92 * G))),
    ]
    for a, b in connectors:
      if b > a:
        micro.append((a, b))

    if len(micro) > budget:
      micro = micro[:budget]
    return micro

  # Run two deterministic micro phases with distinct window families
  # Phase A
  if len(T) <= CAP - CAP_MARGIN:
    room = CAP - len(T)
    if room > 8:
      windows_A = [(0.11, 0.21), (0.34, 0.44), (0.57, 0.67), (0.79, 0.89)]
      microA = build_micro_round(T, room, seed_shift=1, windows=windows_A)
      if microA:
        T.extend(microA[:room])

  if len(T) >= CAP - CAP_MARGIN:
    return T

  # Connector layer B: add a small set of long-range cross-scale ties
  if len(T) <= CAP - CAP_MARGIN:
    T = add_fractional_connectors(T, max_add=3)

  if len(T) >= CAP - CAP_MARGIN:
    return T

  # Phase B (alternate window family)
  if len(T) <= CAP - CAP_MARGIN:
    room = CAP - len(T)
    if room > 8:
      windows_B = [(0.05, 0.15), (0.28, 0.38), (0.60, 0.70), (0.82, 0.92)]
      microB = build_micro_round(T, room, seed_shift=2, windows=windows_B)
      if microB:
        T.extend(microB[:room])

  # Final trim to CAP
  if len(T) > CAP:
    T = T[:CAP]
  return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()