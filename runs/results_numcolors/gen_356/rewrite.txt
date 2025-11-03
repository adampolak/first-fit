# EVOLVE-BLOCK-START

def construct_intervals(seed_count=1):
  """
  Deterministic six-template KT-style spine with two-pass micro phases and long-range connectors.
  Returns a list of open intervals (l, r) in FirstFit presentation order.
  Reference inspiration: Bosek et al., arXiv:1506.00192.
  """

  # Keep total intervals strictly below 10000
  CAP = 9800

  # Fixed, deterministic six-template backbone (rotated per round index)
  template_bank = [
    (2, 6, 10, 14),  # classic KT
    (1, 5, 9, 13),   # left-shifted
    (3, 7, 11, 15),  # right-shifted
    (4, 8, 12, 16),  # stretched-right
    (2, 4, 8, 12),   # compressed inner-left
    (1, 7, 11, 15),  # wide skew
  ]

  # Seed with a single unit interval to avoid early omega inflation
  T = [(0, 1)]

  # Helpers
  def _span_delta(current_T):
    lo = min(l for l, r in current_T)
    hi = max(r for l, r in current_T)
    delta = hi - lo
    if delta <= 0:
      delta = 1
    return lo, hi, delta

  def _append_connectors(S, starts, delta):
    # Classic four connectors that couple blocks without exploding omega
    s0, s1, s2, s3 = starts
    S.append(((s0 - 1) * delta, (s1 - 1) * delta))  # left cap
    S.append(((s2 + 2) * delta, (s3 + 2) * delta))  # right cap
    S.append(((s0 + 2) * delta, (s2 - 1) * delta))  # cross 1
    S.append(((s1 + 2) * delta, (s3 - 1) * delta))  # cross 2

  def _apply_round(current_T, starts, interleave=False, reverse_order=False):
    lo, hi, delta = _span_delta(current_T)

    # Build four translated blocks
    blocks = []
    for s in starts:
      base = s * delta - lo
      block = [(l + base, r + base) for (l, r) in current_T]
      blocks.append(block)

    # Assemble S with optional interleaving and reverse order
    S = []
    if interleave:
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

    # Append connectors
    _append_connectors(S, starts, delta)
    return S

  # Stage 1: six KT rounds with deterministic rotation and interleaving policy
  # Capacity-safe depth for 1-seed, 4-copy + 4-connector recurrence
  depth = 6
  interleave_schedule = [True, False, True, False, True, False]
  reverse_schedule = [False, True, False, True, False, True]

  for ridx in range(depth):
    # Predict next size: sz -> 4*sz + 4
    if 4 * len(T) + 4 > CAP:
      break
    starts = template_bank[ridx % len(template_bank)]
    T = _apply_round(
      T,
      starts,
      interleave=interleave_schedule[ridx % len(interleave_schedule)],
      reverse_order=reverse_schedule[ridx % len(reverse_schedule)]
    )
    if len(T) >= CAP:
      return T[:CAP]

  # Early exit guard
  if len(T) >= CAP - 8:
    return T

  # Stage 1b: deterministic long-range connectors after spine completion
  lo_s, hi_s, delta_s = _span_delta(T)
  G = max(1, hi_s - lo_s)
  def frac_iv(a, b):
    L = lo_s + max(1, int(round(a * G)))
    R = lo_s + max(1, int(round(b * G)))
    if R <= L:
      R = L + 1
    return (L, R)

  # Conservative set of long connectors to reinforce cross-scale coupling
  post_spine_connectors = [
    frac_iv(0.08, 0.60),
    frac_iv(0.25, 0.75),
    frac_iv(0.44, 0.78),
    frac_iv(0.60, 0.92),
  ]
  room = CAP - len(T)
  if room > 0:
    add = post_spine_connectors[:room]
    T.extend(add)

  if len(T) >= CAP - 16:
    return T

  # Two-pass micro phases: CAP-aware, thin seeds, distinct window families
  def _thin_seed(current_T, max_seed):
    n = len(current_T)
    if n <= 0 or max_seed <= 0:
      return []
    stride = max(1, n // max_seed)
    return [current_T[i] for i in range(0, n, stride)][:max_seed]

  def _micro_pass(current_T, budget, pass_id=0):
    if not current_T or budget <= 8:
      return []

    glo = min(l for l, r in current_T)
    ghi = max(r for l, r in current_T)
    G = max(1, ghi - glo)

    # Thin seed, slightly pass-dependent size
    base_seed = max(8, min(40, len(current_T) // (250 if pass_id == 0 else 300)))
    U = _thin_seed(current_T, base_seed)
    if not U:
      return []
    ulo = min(l for l, r in U)

    # Window families (non-overlapping-ish ranges)
    if pass_id == 0:
      windows = [(0.12, 0.22), (0.35, 0.45), (0.58, 0.68), (0.80, 0.90)]
    else:
      windows = [(0.05, 0.15), (0.28, 0.38), (0.60, 0.70), (0.82, 0.92)]

    # Build translated micro-blocks and interleave in a pass-dependent order
    blocks = []
    for (fa, fb) in windows:
      win_lo = glo + int(round(fa * G))
      base = win_lo - ulo
      block = [(l + base, r + base) for (l, r) in U]
      blocks.append(block)

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

    # Fractional-span connectors at the micro scale
    micro_connectors = [
      (glo + int(round(0.08 * G)), glo + int(round(0.30 * G))),  # left cap
      (glo + int(round(0.26 * G)), glo + int(round(0.56 * G))),  # cross 1
      (glo + int(round(0.44 * G)), glo + int(round(0.78 * G))),  # cross 2
      (glo + int(round(0.60 * G)), glo + int(round(0.92 * G))),  # right cap
    ]
    for a, b in micro_connectors:
      if b > a:
        micro.append((a, b))

    # Trim to available budget
    if len(micro) > budget:
      micro = micro[:budget]
    return micro

  # Split remaining capacity into two passes
  remaining = CAP - len(T)
  if remaining > 8:
    budget_A = remaining // 2
    A = _micro_pass(T, budget_A, pass_id=0)
    if A:
      T.extend(A)

  remaining = CAP - len(T)
  if remaining > 8:
    budget_B = remaining
    B = _micro_pass(T, budget_B, pass_id=1)
    if B:
      if len(B) > remaining:
        B = B[:remaining]
      T.extend(B)

  # Final trim to CAP
  if len(T) > CAP:
    T = T[:CAP]

  return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()