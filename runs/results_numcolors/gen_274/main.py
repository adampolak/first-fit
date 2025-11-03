# EVOLVE-BLOCK-START

def construct_intervals(seed_count=1):
  """
  Construct a sequence of intervals on the real line, in the order presented to FirstFit,
  to maximize FF colors divided by the clique number.
  Returns:
    intervals: list of (l, r) tuples (open intervals).
  """

  # Capacity guard tuned to return after six spine rounds (skipping micro-phase) and < 10000
  CAP = 9800

  # Stage 1: Deterministic KT-style spine (no interleaving), 6 rounds, classic connectors.
  spine_starts = (2, 6, 10, 14)

  # Seed: keep a single unit interval; multi-seed tends to inflate omega too early.
  T = [(0, 1)] if seed_count <= 1 else [(0, 1)]

  def apply_spine_round(current_T, starts, do_interleave=False, reverse_order=False):
    # Compute span and delta
    lo = min(l for l, r in current_T)
    hi = max(r for l, r in current_T)
    delta = hi - lo
    if delta <= 0:
      delta = 1

    # Build translated blocks
    blocks = []
    for s in starts:
      base = s * delta - lo
      block = [(l + base, r + base) for (l, r) in current_T]
      blocks.append(block)

    # Assemble S with optional interleaving and order control
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

    # Classic connectors (Figure 4 style)
    s0, s1, s2, s3 = starts
    connectors = [
      ((s0 - 1) * delta, (s1 - 1) * delta),  # left cap
      ((s2 + 2) * delta, (s3 + 2) * delta),  # right cap
      ((s0 + 2) * delta, (s2 - 1) * delta),  # cross 1
      ((s1 + 2) * delta, (s3 - 1) * delta),  # cross 2
    ]
    S.extend(connectors)
    return S

  # Perform exactly six spine rounds, rotating start templates to diversify structure.
  template_bank = [
    (2, 6, 10, 14),  # classic KT
    (1, 5, 9, 13),   # left-shifted
    (3, 7, 11, 15),  # right-shifted
    (4, 8, 12, 16),  # alternate
  ]
  bank_size = len(template_bank)
  for ridx in range(6):
    # Rotate starts per round to break symmetry
    starts = template_bank[ridx % bank_size]
    # Predict next size: sz -> 4*sz + 4
    nxt_size = 4 * len(T) + 4
    if nxt_size > CAP:
      break
    # Parity-based interleaving: even rounds interleave; odd rounds reverse order
    do_interleave = (ridx % 2 == 0)
    reverse_order = (ridx % 2 == 1)
    T = apply_spine_round(T, starts, do_interleave=do_interleave, reverse_order=reverse_order)

  # If we are already near the cap, return the strong baseline.
  # Inject a tiny set of long caps near the tail before micro-phases to amplify FF mixing.
  lo = min(l for l, r in T)
  hi = max(r for l, r in T)
  span = max(1, hi - lo)
  def _cap_at(a_frac, b_frac):
    L = lo + max(1, int(round(a_frac * span)))
    R = lo + max(1, int(round(b_frac * span)))
    if R <= L:
      R = L + 1
    return (L, R)
  tail_caps = [_cap_at(0.08, 0.60), _cap_at(0.25, 0.75), _cap_at(0.75, 0.92)]
  for c in tail_caps:
    if len(T) >= CAP:
      break
    T.append(c)
  if len(T) >= CAP - 8:
    return T

  # Stage 2: Fractional-window micro extension rounds (thin sampling; interleaved).
  # Goals:
  #  - Raise FF pressure via cross-scale coupling and interleaving parity,
  #  - Keep omega in check by using thin seeds and fractional-span connectors,
  #  - Respect strict capacity guard.

  def thin_seed(current_T, max_seed):
    """Take a thin, evenly spaced sample of current_T of size <= max_seed."""
    n = len(current_T)
    if n == 0 or max_seed <= 0:
      return []
    step = max(1, n // max_seed)
    U = current_T[::step][:max_seed]
    return U

  def micro_round(current_T, round_id, budget):
    if budget <= 0 or not current_T:
      return []

    glo = min(l for l, r in current_T)
    ghi = max(r for l, r in current_T)
    G = max(1, ghi - glo)

    # Thin seed: bounded and deterministic size to respect budget
    per_block_target = max(16, min(40, budget // 10))
    U = thin_seed(current_T, per_block_target)
    if not U:
      return []

    ulo = min(l for l, r in U)

    # Four fractional windows with a small parity-based shift
    shift = (round_id % 3) * 0.02
    window_fracs = [
      (0.12 + shift, 0.22 + shift),
      (0.35 + shift, 0.45 + shift),
      (0.58 + shift, 0.68 + shift),
      (0.80 + shift, 0.90 + shift),
    ]
    window_fracs = [
      (max(0.05, min(0.90, a)), max(0.10, min(0.95, b)))
      for (a, b) in window_fracs
    ]

    # Build translated micro-blocks aligned to these windows
    blocks = []
    for idx, (fa, fb) in enumerate(window_fracs):
      win_lo = glo + int(round(fa * G))
      base = win_lo - ulo
      block = [(l + base, r + base) for (l, r) in U]
      # Alternate internal reversal by window index and round parity
      if ((idx + round_id) % 2) == 1:
        block = list(reversed(block))
      blocks.append(block)

    # Interleave micro-blocks (reverse order on odd rounds)
    micro = []
    maxlen = max(len(b) for b in blocks)
    order = list(range(len(blocks)))
    if round_id % 2 == 1:
      order.reverse()
    for i in range(maxlen):
      for j in order:
        blk = blocks[j]
        if i < len(blk):
          micro.append(blk[i])

    # Fractional-span connectors at the micro scale
    micro_connectors = [
      (glo + int(round(0.08 * G)), glo + int(round(0.30 * G))),  # left cap
      (glo + int(round(0.60 * G)), glo + int(round(0.92 * G))),  # right cap
      (glo + int(round(0.26 * G)), glo + int(round(0.56 * G))),  # cross 1
      (glo + int(round(0.44 * G)), glo + int(round(0.78 * G))),  # cross 2
    ]
    for a, b in micro_connectors:
      if b > a:
        micro.append((a, b))

    # Enforce budget
    if len(micro) > budget:
      micro = micro[:budget]
    return micro

  # Execute two micro rounds with strong guards
  remaining = CAP - len(T)
  mr1 = micro_round(T, round_id=0, budget=max(0, remaining // 2))
  if mr1:
    room = CAP - len(T)
    if len(mr1) > room:
      mr1 = mr1[:room]
    T.extend(mr1)

  remaining = CAP - len(T)
  mr2 = micro_round(T, round_id=1, budget=max(0, remaining))
  if mr2:
    room = CAP - len(T)
    if len(mr2) > room:
      mr2 = mr2[:room]
    T.extend(mr2)

  return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()