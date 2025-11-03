# EVOLVE-BLOCK-START

def construct_intervals(enable_alt_microphase=True):
  """
  Construct a sequence of intervals on the real line, in the order presented to FirstFit,
  to maximize FF colors divided by the clique number (target omega <= 10),
  using a two-pass CAP-aware design:
    - Pass 1: deterministic KT-like spine with a bounded shrinker
    - Pass 2: two micro phases (window families) with hard budgeting and pin staggering

  Args:
    enable_alt_microphase (bool): enable the second micro-phase with alternate windows.

  Returns:
    intervals: list of tuples (l, r) representing open intervals.
  """

  # Hard capacity to keep total intervals < 10000
  CAP = 9800

  # Pass-1: deterministic KT-like spine templates (six-template rotation)
  spine_templates = [
    (2, 6, 10, 14),  # classic KT
    (1, 5, 9, 13),   # left-shifted
    (3, 7, 11, 15),  # right-shifted
    (4, 8, 12, 16),  # stretched-right
    (5, 9, 13, 17),  # extra right-shift
    (6, 10, 14, 18), # extra classic-shift
  ]

  # Seed: a single unit interval to keep omega low at the start
  T = [(0, 1)]

  # Helpers
  def span_of(seq):
    lo = min(l for l, r in seq)
    hi = max(r for l, r in seq)
    d = hi - lo
    if d <= 0:
      d = 1
    return lo, hi, d

  def kt_round(current_T, starts, interleave=False, reverse=False):
    lo, hi, delta = span_of(current_T)

    # Build translated blocks
    blocks = []
    for s in starts:
      base = s * delta - lo
      blocks.append([(l + base, r + base) for (l, r) in current_T])

    # Assemble S with policy
    S = []
    if interleave:
      order = list(range(4))
      if reverse:
        order.reverse()
      maxlen = max(len(b) for b in blocks)
      for i in range(maxlen):
        for idx in order:
          blk = blocks[idx]
          if i < len(blk):
            S.append(blk[i])
    else:
      if reverse:
        blocks = list(reversed(blocks))
      for blk in blocks:
        S.extend(blk)

    # Classic four connectors
    s0, s1, s2, s3 = starts
    S.append(((s0 - 1) * delta, (s1 - 1) * delta))
    S.append(((s2 + 2) * delta, (s3 + 2) * delta))
    S.append(((s0 + 2) * delta, (s2 - 1) * delta))
    S.append(((s1 + 2) * delta, (s3 - 1) * delta))
    return S

  # Per-round bounded shrinker: carve a few thin pins from a small uniform seed
  def bounded_shrinker(seq, cap_left):
    if cap_left <= 0 or not seq:
      return []
    lo, hi, d = span_of(seq)
    # Thin uniform seed: up to 12 samples or 1% of size
    max_samples = max(6, min(12, max(1, len(seq) // 100)))
    stride = max(1, len(seq) // max_samples)
    seed = [seq[i] for i in range(0, len(seq), stride)][:max_samples]
    if not seed:
      return []
    # Pin length ~ d/2000, at least 1
    eps = max(1, d // 2000)
    out = []
    # Deterministic 3-class staggering to reduce local overlap
    for idx, (l, r) in enumerate(seed):
      if r - l <= 2 * eps:
        continue
      mid = (l + r) // 2
      shift = (idx % 3) * (eps // 2 + 1)
      L = mid - eps + shift
      R = L + eps
      if R > L:
        out.append((L, R))
      if len(out) >= cap_left:
        break
    return out

  # Execute spine with deterministic parity interleaving and shrinker
  for ridx in range(6):  # six rounds; each round ~ x4 + 4
    if 4 * len(T) + 4 > CAP:
      break
    starts = spine_templates[ridx % len(spine_templates)]
    interleave = (ridx % 2 == 0)   # interleave on even rounds
    reverse = (ridx % 2 == 1)      # reverse block order on odd rounds
    T = kt_round(T, starts, interleave, reverse)

    # Insert a few thin pins (bounded) to increase FF pressure safely
    room = CAP - len(T)
    if room > 0:
      pins = bounded_shrinker(T, min(room, 9))  # at most 9 pins per round
      if pins:
        # Insert near the end to affect late FF decisions
        T.extend(pins)

  if len(T) >= CAP - 24:
    return T

  # Dedicated long-range connector augmentation pass (4–6 connectors)
  lo_s, hi_s, d_s = span_of(T)
  def frac(a):
    return lo_s + max(1, int(round(a * d_s)))
  connectors = [
    (frac(0.06), frac(0.28)),
    (frac(0.18), frac(0.46)),
    (frac(0.34), frac(0.66)),
    (frac(0.54), frac(0.84)),
    (frac(0.72), frac(0.94)),
    (frac(0.22), frac(0.78)),
  ]
  for a, b in connectors:
    if len(T) >= CAP:
      break
    if b > a:
      T.append((a, b))

  if len(T) >= CAP - 24:
    return T

  # Pass-2: Two micro phases with hard CAP budgeting and omega-conscious pinning

  # Budgeting: reserve a fixed portion for phase-1 windows and the rest for phase-2
  remaining = CAP - len(T)
  if remaining <= 0:
    return T

  # Deterministic split: 65% to phase-1, 35% to phase-2 (if enabled)
  p1_budget = int(0.65 * remaining)
  p2_budget = remaining - p1_budget if enable_alt_microphase else 0

  # Phase-1 window family (central windows)
  phase1_windows = [(0.12, 0.22), (0.35, 0.45), (0.58, 0.68), (0.80, 0.90)]
  # Phase-2 window family (alternate windows)
  phase2_windows = [(0.05, 0.15), (0.28, 0.38), (0.60, 0.70), (0.82, 0.92)]

  def micro_phase(current_T, budget, windows, with_long_cross=False):
    """Build a micro phase from a thin seed, translated into given windows,
       transformed into thin pins with staggering and guarded by a local cap."""
    if budget <= 0 or not current_T:
      return []

    glo = min(l for l, r in current_T)
    ghi = max(r for l, r in current_T)
    G = max(1, ghi - glo)

    # Thin, evenly spaced seed
    seed_sz = max(8, min(36, len(current_T) // 280))
    stride = max(1, len(current_T) // max(1, seed_sz))
    U = [current_T[i] for i in range(0, len(current_T), stride)][:seed_sz]
    if not U:
      return []

    ulo = min(l for l, r in U)
    # Pin length at micro scale
    eps = max(1, G // 4000)

    # Local concurrency governor per window: cap pins per window to keep omega small
    per_window_cap = max(8, min(20, budget // max(1, len(windows) * 4)))

    blocks = []
    for wid, (fa, fb) in enumerate(windows):
      win_lo = glo + int(round(fa * G))
      base = win_lo - ulo
      block = []
      for idx, (l, r) in enumerate(U):
        # Convert to a pin with deterministic stagger inside the window
        mid = (l + r) // 2
        shift = (idx % 4) * (eps // 2 + 1) + wid  # window-coupled stagger
        L = mid + base - (eps // 2) + shift
        R = L + eps
        if R > L:
          block.append((L, R))
        if len(block) >= per_window_cap:  # governor
          break
      blocks.append(block)

    # Interleave blocks deterministically
    M = []
    maxlen = max((len(b) for b in blocks), default=0)
    for i in range(maxlen):
      for blk in blocks:
        if i < len(blk):
          M.append(blk[i])

    # Deterministic micro-connectors across windows
    conn = [
      (glo + int(round(0.08 * G)), glo + int(round(0.30 * G))),
      (glo + int(round(0.26 * G)), glo + int(round(0.56 * G))),
      (glo + int(round(0.44 * G)), glo + int(round(0.78 * G))),
      (glo + int(round(0.60 * G)), glo + int(round(0.92 * G))),
    ]
    for a, b in conn:
      if b > a:
        M.append((a, b))

    # Optional long cross to tie distant colors
    if with_long_cross:
      a = glo + int(round(0.18 * G))
      b = glo + int(round(0.84 * G))
      if b > a:
        M.append((a, b))

    # Trim to budget
    if len(M) > budget:
      M = M[:budget]
    return M

  # Execute phase-1
  if p1_budget > 0:
    phase1 = micro_phase(T, p1_budget, phase1_windows, with_long_cross=False)
    if phase1:
      T.extend(phase1)

  # Recompute remaining budget
  remaining = CAP - len(T)
  if remaining <= 0:
    return T

  # Execute phase-2 (optional)
  if enable_alt_microphase and p2_budget > 0:
    # Recompute per current T to anchor to new span
    phase2 = micro_phase(T, min(p2_budget, CAP - len(T)), phase2_windows, with_long_cross=True)
    if phase2:
      T.extend(phase2)

  # Final capacity trim
  if len(T) > CAP:
    T = T[:CAP]

  return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()