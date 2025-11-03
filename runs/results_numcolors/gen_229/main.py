# EVOLVE-BLOCK-START

def construct_intervals(seed_count=1):
  """
  Construct a sequence of intervals on the real line, in the order presented to FirstFit,
  to maximize FirstFit colors divided by the clique number (omega).

  Structural pipeline:
    1) KT-style rotating-template backbone (six rounds, parity interleaving)
    2) Sparse long caps layer (tail insert)
    3) Micro-phase A (primary fractional windows)
    4) Micro-phase B (alternate fractional windows)
    5) Cross4 connector layer (long-range couplers, gated)

  Inputs:
    seed_count (int): kept for API compatibility; multi-seeding is guarded.

  Output:
    list of (l, r) integer tuples (open intervals), in presentation order.
  """

  # ------------------------------
  # Global guards and parameters
  # ------------------------------
  CAP = 9800  # hard capacity guard
  BASE_SEED = 911382323

  # Templates for the KT spine
  TEMPLATE_BANK = [
    (2, 6, 10, 14),  # classic KT
    (1, 5, 9, 13),   # left-shifted
    (3, 7, 11, 15),  # right-shifted
    (4, 8, 12, 16),  # stretched-right
  ]

  # Window families for micro-phases
  WINDOWS_PRIMARY = [(0.12, 0.22), (0.35, 0.45), (0.58, 0.68), (0.80, 0.90)]
  WINDOWS_ALT     = [(0.05, 0.15), (0.28, 0.38), (0.60, 0.70), (0.82, 0.92)]

  # ------------------------------
  # Utility helpers
  # ------------------------------
  def _span(T):
    lo = min(l for l, _ in T)
    hi = max(r for _, r in T)
    G = hi - lo
    if G <= 0:
      G = 1
    return lo, hi, G

  def _predict_next_size(sz):
    # KT recurrence: n -> 4*n + 4 (four copies + four classic connectors)
    return 4 * sz + 4

  def _dhash(*xs):
    # Small deterministic hash, stable across runs
    h = BASE_SEED
    for v in xs:
      h = ((h ^ (v + 0x9e3779b97f4a7c15)) * 0xbf58476d1ce4e5b9) & 0xFFFFFFFFFFFFFFFF
      h ^= (h >> 27)
    return h

  def _thin_seed(T, max_seed):
    n = len(T)
    if n == 0 or max_seed <= 0:
      return []
    step = max(1, n // max_seed)
    return T[::step][:max_seed]

  # ------------------------------
  # Spine builder (KT-style rounds)
  # ------------------------------
  def _apply_spine_round(T, starts, interleave=False, reverse_order=False):
    lo, hi, G = _span(T)
    blocks = []
    for s in starts:
      base = s * G - lo
      blocks.append([(l + base, r + base) for (l, r) in T])

    S = []
    if interleave:
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

    s0, s1, s2, s3 = starts
    S.append(((s0 - 1) * G, (s1 - 1) * G))  # left cap
    S.append(((s2 + 2) * G, (s3 + 2) * G))  # right cap
    S.append(((s0 + 2) * G, (s2 - 1) * G))  # cross 1
    S.append(((s1 + 2) * G, (s3 - 1) * G))  # cross 2
    return S

  def build_spine(rounds=6, rotate=True):
    # Seed: single unit interval to keep early omega low
    T = [(0, 1)] if seed_count <= 1 else [(0, 1)]
    for ridx in range(max(0, int(rounds))):
      nxt = _predict_next_size(len(T))
      if nxt > CAP:
        break
      starts = TEMPLATE_BANK[ridx % len(TEMPLATE_BANK)] if rotate else TEMPLATE_BANK[0]
      inter = (ridx % 2 == 0)
      rev = (ridx % 2 == 1)
      T = _apply_spine_round(T, starts, interleave=inter, reverse_order=rev)
      if len(T) >= CAP:
        return T[:CAP]
    return T

  # ------------------------------
  # Caps layer (sparse long caps)
  # ------------------------------
  def inject_caps(T):
    if len(T) >= CAP - 8:
      return T
    lo, hi, G = _span(T)
    def cap_at(a, b):
      L = lo + max(1, int(round(a * G)))
      R = lo + max(1, int(round(b * G)))
      if R <= L:
        R = L + 1
      return (L, R)
    caps = [cap_at(0.08, 0.60), cap_at(0.25, 0.75), cap_at(0.75, 0.92)]
    for c in caps:
      if len(T) >= CAP:
        break
      T.append(c)
    return T

  # ------------------------------
  # Micro-phase builders
  # ------------------------------
  def _build_micro_windows(T, budget, windows, iter_id=0, jitter=False):
    if not T or budget <= 8:
      return []
    glo, ghi, G = _span(T)
    # Thin, evenly spaced seed; bounded size for safety
    seed_sz = max(8, min(40, len(T) // 250))
    stride = max(1, len(T) // max(1, seed_sz))
    U = [T[i] for i in range(0, len(T), stride)][:seed_sz]
    if not U:
      return []
    ulo = min(l for l, _ in U)

    # Optional tiny deterministic jitter across windows
    j = 0.0
    if jitter:
      h = _dhash(iter_id, len(T))
      j = ((h % 3) - 1) * 0.005  # in {-0.005, 0, 0.005}

    blocks = []
    for (fa, fb) in windows:
      a = max(0.05, min(0.90, fa + j))
      b = max(a + 0.02, min(0.95, fb + j))
      win_lo = glo + int(round(a * G))
      base = win_lo - ulo
      block = [(l + base, r + base) for (l, r) in U]
      tag = ((int(round(a * 100)) // 5) + iter_id) % 2
      if tag == 1:
        block = list(reversed(block))
      blocks.append(block)

    # Interleave with parity-based order
    order = list(range(len(blocks)))
    if iter_id % 2 == 1:
      order.reverse()
    micro = []
    maxlen = max(len(b) for b in blocks)
    for i in range(maxlen):
      for idx in order:
        blk = blocks[idx]
        if i < len(blk):
          micro.append(blk[i])

    # Fractional connectors (KT-analog at global scale)
    connectors = [
      (glo + int(round(0.08 * G)), glo + int(round(0.30 * G))),
      (glo + int(round(0.60 * G)), glo + int(round(0.92 * G))),
      (glo + int(round(0.26 * G)), glo + int(round(0.56 * G))),
      (glo + int(round(0.44 * G)), glo + int(round(0.78 * G))),
    ]
    for a, b in connectors:
      if b > a:
        micro.append((a, b))

    # Budget trim
    if len(micro) > budget:
      micro = micro[:budget]
    return micro

  # ------------------------------
  # Cross4 connectors (long-range ties)
  # ------------------------------
  def cross4_layer(T, max_add=12):
    if not T or max_add <= 0:
      return T
    glo, ghi, G = _span(T)
    # Deterministic placement, conservative count and CAP-gated
    cross = [
      (glo + int(round(0.10 * G)), glo + int(round(0.85 * G))),
      (glo + int(round(0.18 * G)), glo + int(round(0.72 * G))),
      (glo + int(round(0.32 * G)), glo + int(round(0.94 * G))),
    ]
    for iv in cross[:max_add]:
      if len(T) >= CAP:
        break
      if iv[1] > iv[0]:
        T.append(iv)
    return T

  # ------------------------------
  # Orchestrate pipeline
  # ------------------------------
  T = build_spine(rounds=6, rotate=True)

  if len(T) < CAP - 8:
    T = inject_caps(T)

  # Micro-phase A (primary windows), CAP-aware
  room = CAP - len(T)
  if room > 8:
    microA = _build_micro_windows(T, budget=min(room, max(12, room // 2)),
                                  windows=WINDOWS_PRIMARY, iter_id=0, jitter=False)
    if microA:
      take = min(len(microA), CAP - len(T))
      T.extend(microA[:take])

  # Micro-phase B (alternate windows), CAP-aware
  room = CAP - len(T)
  if room > 8:
    microB = _build_micro_windows(T, budget=min(room, 128),
                                  windows=WINDOWS_ALT, iter_id=1, jitter=False)
    if microB:
      take = min(len(microB), CAP - len(T))
      T.extend(microB[:take])

  # Cross4 long-range connector layer (final small push), CAP-aware
  room = CAP - len(T)
  if room > 0:
    T = cross4_layer(T, max_add=min(8, room))

  # Final capacity trim
  if len(T) > CAP:
    T = T[:CAP]

  return T


def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()