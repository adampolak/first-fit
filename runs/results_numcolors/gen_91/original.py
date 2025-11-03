# EVOLVE-BLOCK-START

def construct_intervals(rounds=6,
                        rotate_starts=True,
                        reverse_block_parity=True,
                        interleave_blocks=True,
                        phase2_iters=1):
  """
  Deterministic KT-style spine with a safeguarded two-phase scaffold.

  Parameters:
    rounds (int): main expansion depth; near 6 yields ~9556 intervals for a single-seed KT spine.
    rotate_starts (bool): available but disabled in the safeguarded spine mode.
    reverse_block_parity (bool): available but disabled in the safeguarded spine mode.
    interleave_blocks (bool): available but disabled in the safeguarded spine mode.
    phase2_iters (int): requested micro iterations; actual application is strictly capacity- and safety-gated.

  Returns:
    intervals: list of (l, r) integer tuples, open intervals, in FF presentation order.
  """

  # Hard capacity guard to keep the total count < 10000
  CAP = 9800

  # Stable KT-spine start pattern (proven strong in prior runs)
  spine_starts = [2, 6, 10, 14]

  # Optional template bank for exploratory use in micro/secondary phase
  # (kept dormant by default to preserve the strong baseline).
  template_bank = [
    [2, 6, 10, 14],  # A: classic KT
    [1, 5, 9, 13],   # B: left-shifted
    [3, 7, 11, 15],  # C: right-shifted
    [2, 4, 8, 12],   # D: compressed left pair
    [3, 5, 9, 13],   # E: gentle left pack
    [1, 7, 11, 15],  # F: wide skew
    [2, 8, 10, 12],  # G: inner symmetric
    [4, 6, 8, 10],   # H: tight middle
  ]

  # Seed with one unit interval to allow six KT rounds within CAP.
  T = [(0, 1)]

  # Predictive size accounting to cap the number of full rounds.
  # KT growth per round: size -> 4*size + 4
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

  # Stage 1: KT spine, un-interleaved, classic connectors, capacity-safe depth.
  depth, size_est = max_rounds_within_cap(len(T), rounds)

  def apply_round(current_T, starts):
    # Compute span and delta
    lo = min(l for l, r in current_T)
    hi = max(r for l, r in current_T)
    delta = hi - lo
    if delta <= 0:
      delta = 1

    # Build four translated blocks sequentially (no interleaving)
    S = []
    for s in starts:
      base = s * delta - lo
      # No parity reversal here: spine mode
      S.extend((l + base, r + base) for (l, r) in current_T)

    # Classic connectors that preserve the desired FF pressure without
    # blowing up omega when used with spine_starts:
    s0, s1, s2, s3 = starts
    connectors = [
      ((s0 - 1) * delta, (s1 - 1) * delta),
      ((s2 + 2) * delta, (s3 + 2) * delta),
      ((s0 + 2) * delta, (s2 - 1) * delta),
      ((s1 + 2) * delta, (s3 - 1) * delta),
    ]
    S.extend(connectors)
    return S

  for _ in range(depth):
    T = apply_round(T, spine_starts)

  # If we already reached or are near the cap, skip phase 2 safely.
  if len(T) >= CAP - 16:
    return T

  # Stage 2 (safeguarded): attempt up to phase2_iters micro steps using a smaller delta2.
  # To avoid harming the strong KT baseline, we:
  #  - require enough residual capacity,
  #  - use at most one light-touch micro step by default,
  #  - place micro additions as sparse long caps that do not densify the core "spine."
  def micro_caps(current_T, intensity=1):
    if not current_T:
      return []
    lo = min(l for l, r in current_T)
    hi = max(r for l, r in current_T)
    delta = max(1, hi - lo)
    d2 = max(1, delta // 3)

    # Two very sparse long caps placed away from the densest overlap
    caps = [
      (lo + 1 * d2, lo + 5 * d2),
      (hi - 6 * d2, hi - 2 * d2),
    ]
    # intensity>1 adds a symmetric mid cap pair (still light)
    if intensity > 1:
      mid = (lo + hi) // 2
      caps += [
        (mid - 2 * d2, mid + 2 * d2),
      ]
    return caps

  # Compute how many safe micro intervals we can add
  # Respect requested phase2_iters, but safeguard with CAP and keep additions tiny.
  micro_steps = max(0, int(phase2_iters))
  # Each step adds at most 2–3 caps; plan conservatively.
  add_per_step = 2
  max_steps_by_cap = max(0, (CAP - len(T)) // add_per_step)

  steps = min(micro_steps, max_steps_by_cap, 1)  # at most 1 guarded micro pass by default
  for k in range(steps):
    # keep micro intensity very low to avoid increasing omega
    caps = micro_caps(T, intensity=1 if k == 0 else 1)
    # capacity guard
    room = CAP - len(T)
    if room <= 0:
      break
    if len(caps) > room:
      caps = caps[:room]
    T.extend(caps)

  return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()