# EVOLVE-BLOCK-START

def construct_intervals(seed_count=1):
  """
  Construct a sequence of intervals on the real line, in the order presented to FirstFit,
  to maximize FF colors divided by the clique number (omega).

  Interface preserved: construct_intervals(seed_count=1) -> list[(l, r)]
  """

  # Hard capacity guard to keep total intervals < 10000
  CAP = 9800

  # Deterministic base schedule (no PRNG; fully reproducible)
  BASE_SEED = 1729  # reserved for potential hashing; schedule here is purely index-driven

  # Primary rotating four-start templates for the KT spine
  template_bank = [
    (2, 6, 10, 14),  # classic KT
    (1, 5, 9, 13),   # left-shifted
    (3, 7, 11, 15),  # right-shifted
    (4, 8, 12, 16),  # stretched-right
  ]

  # Auxiliary template for a light density boost (kept modest)
  density_template = (2, 4, 8, 12)

  # Seed with a single unit interval; multi-seed tends to inflate omega too early
  T = [(0, 1)]

  # ---------- Helpers ----------

  def _span_delta(current_T):
    lo = min(l for l, r in current_T)
    hi = max(r for l, r in current_T)
    delta = hi - lo
    if delta <= 0:
      delta = 1
    return lo, hi, delta

  def _append_connectors(S, starts, delta):
    # Classic four connectors; preserves strong FF pressure while keeping omega modest.
    s0, s1, s2, s3 = starts
    S.append(((s0 - 1) * delta, (s1 - 1) * delta))  # left cap
    S.append(((s2 + 2) * delta, (s3 + 2) * delta))  # right cap
    S.append(((s0 + 2) * delta, (s2 - 1) * delta))  # cross 1
    S.append(((s1 + 2) * delta, (s3 - 1) * delta))  # cross 2

  def _apply_round(current_T, starts, do_interleave=False, reverse_order=False):
    lo, hi, delta = _span_delta(current_T)

    # Build four translated blocks
    blocks = []
    for s in starts:
      base = s * delta - lo
      block = [(l + base, r + base) for (l, r) in current_T]
      blocks.append(block)

    # Assemble with optional interleaving and reverse order to mix colors
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

    # Append classic connectors
    _append_connectors(S, starts, delta)
    return S

  # Deterministic, CAP-aware thin seed picker
  def _thin_seed(current_T, max_seed):
    n = len(current_T)
    if n <= 0 or max_seed <= 0:
      return []
    step = max(1, n // max_seed)
    return current_T[::step][:max_seed]

  # Micro-round builder with custom window set
  def _build_micro_round(current_T, budget, window_fracs):
    if not current_T or budget <= 8:
      return []

    glo = min(l for l, r in current_T)
    ghi = max(r for l, r in current_T)
    G = max(1, ghi - glo)

    # Thin, evenly spaced seed (kept small to not inflate omega)
    seed_sz = max(8, min(36, len(current_T) // 270))
    U = _thin_seed(current_T, seed_sz)
    if not U:
      return []

    ulo = min(l for l, r in U)

    # Build translated micro-blocks aligned to the provided windows
    blocks = []
    for (fa, fb) in window_fracs:
      # Clamp windows strictly inside (0.05, 0.95) to avoid boundary pileups
      a = max(0.05, min(0.90, fa))
      b = max(0.10, min(0.95, fb))
      win_lo = glo + int(round(a * G))
      base = win_lo - ulo
      block = [(l + base, r + base) for (l, r) in U]
      # Internal reversal by parity of window index improves mixing
      if int(round(a * 100)) // 5 % 2 == 1:
        block = list(reversed(block))
      blocks.append(block)

    # Interleave micro-blocks forward
    micro = []
    maxlen = max(len(b) for b in blocks)
    for i in range(maxlen):
      for blk in blocks:
        if i < len(blk):
          micro.append(blk[i])

    # Fractional-span connectors at the micro scale (keep modest, deterministic)
    micro_connectors = [
      (glo + int(round(0.08 * G)), glo + int(round(0.30 * G))),  # left cap
      (glo + int(round(0.60 * G)), glo + int(round(0.92 * G))),  # right cap
      (glo + int(round(0.26 * G)), glo + int(round(0.56 * G))),  # cross 1
      (glo + int(round(0.44 * G)), glo + int(round(0.78 * G))),  # cross 2
    ]
    for a, b in micro_connectors:
      if b > a:
        micro.append((a, b))

    # Trim to available budget
    if len(micro) > budget:
      micro = micro[:budget]
    return micro

  # Lightweight density boost: echo a thin seed across an auxiliary template at spine scale
  def _density_boost(current_T, starts):
    if not current_T:
      return []
    lo, hi, delta = _span_delta(current_T)
    seed = _thin_seed(current_T, max_seed=24)
    if not seed:
      return []
    ulo = min(l for l, r in seed)
    blocks = []
    # Use a deterministic auxiliary template at the same scale (small)
    aux = density_template
    for s in aux:
      base = s * delta - ulo
      block = [(l + base, r + base) for (l, r) in seed]
      blocks.append(block)
    # Interleave to mix
    boost = []
    maxlen = max(len(b) for b in blocks)
    for i in range(maxlen):
      for blk in blocks:
        if i < len(blk):
          boost.append(blk[i])
    # Append a pair of modest connectors
    s0, s1, s2, s3 = aux
    boost.append(((s0 + 1) * delta, (s2 - 1) * delta))
    boost.append(((s1 + 1) * delta, (s3 - 1) * delta))
    return boost

  # ---------- Stage 1: KT spine, deterministic 4-cycle interleaving ----------

  # Next-size predictor: n -> 4n + 4
  def _next_size(sz):
    return 4 * sz + 4

  # Deterministic per-round interleaving schedule (4-cycle)
  def _round_policy(idx):
    # Cycle: [(interleave, reverse_order)]
    cycle = [(True, False), (False, True), (True, False), (False, False)]
    return cycle[idx % len(cycle)]

  for ridx in range(6):
    if _next_size(len(T)) > CAP:
      break
    starts = template_bank[ridx % len(template_bank)]
    do_interleave, reverse_order = _round_policy(ridx)
    T = _apply_round(T, starts, do_interleave=do_interleave, reverse_order=reverse_order)
    if len(T) >= CAP:
      T = T[:CAP]
      return T

  # Early exit if near capacity
  if len(T) >= CAP - 8:
    return T

  # ---------- Lightweight density boost tied to the spine ----------
  # CAP-aware: allow only a small number to avoid omega inflation
  room = CAP - len(T)
  if room > 16:
    starts_hint = template_bank[(BASE_SEED + len(T)) % len(template_bank)]
    boost = _density_boost(T, starts_hint)
    if boost:
      if len(boost) > room:
        boost = boost[:room]
      T.extend(boost)

  # If still near capacity, return
  if len(T) >= CAP - 8:
    return T

  # ---------- Stage 2: Dual micro-phases with distinct windows ----------

  # Micro-phase A: baseline windows (inside span)
  windows_A = [(0.12, 0.22), (0.35, 0.45), (0.58, 0.68), (0.80, 0.90)]
  room = CAP - len(T)
  if room > 8:
    microA = _build_micro_round(T, room, windows_A)
    if microA:
      if len(microA) > room:
        microA = microA[:room]
      T.extend(microA)

  # Micro-phase B: recommended second window set to couple missed interactions
  windows_B = [(0.16, 0.26), (0.40, 0.50), (0.62, 0.72), (0.84, 0.94)]
  room = CAP - len(T)
  if room > 8:
    microB = _build_micro_round(T, room, windows_B)
    if microB:
      if len(microB) > room:
        microB = microB[:room]
      T.extend(microB)

  # ---------- Deterministic long-range cross-scale connectors ----------
  room = CAP - len(T)
  if room > 0:
    lo2, hi2, delta2 = _span_delta(T)
    def _frac_iv(a_frac, b_frac):
      L = lo2 + max(1, int(round(a_frac * delta2)))
      R = lo2 + max(1, int(round(b_frac * delta2)))
      if R <= L:
        R = L + 1
      return (L, R)
    connectors2 = [
      _frac_iv(0.14, 0.52),
      _frac_iv(0.32, 0.88),
      _frac_iv(0.05, 0.47),
      _frac_iv(0.62, 0.95),
      _frac_iv(0.20, 0.41),
      _frac_iv(0.71, 0.89),
    ]
    if room < len(connectors2):
      connectors2 = connectors2[:room]
    T.extend(connectors2)

  # Final capacity trim
  if len(T) > CAP:
    T = T[:CAP]

  return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()