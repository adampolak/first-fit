# EVOLVE-BLOCK-START

def construct_intervals(rounds=3, seed_lo=0.0, seed_scale=1.0):
  """
  Multiscale, deterministic two-phase construction to pressure FirstFit
  while keeping the offline optimum (clique number) small (target <= 10).

  Interface preserved:
    - rounds: unused as a hard cap; we adaptively choose main rounds to fit budget.
    - seed_lo, seed_scale: initial seed placement and scale.

  Returns:
    A list of (l, r) real pairs representing open intervals in presentation order.
  """

  # Hard cap for practical reasons (evaluator budget)
  MAX_TOTAL = 9800

  # Main-phase parameters
  # We choose 6 rounds as default (empirically strong FF usage under budget with 4-connectors).
  MAIN_ROUNDS = 6

  # Secondary delta2 micro-phase parameters
  MICRO_ROUNDS = 2   # try two micro rounds if budget allows
  MICRO_BLOCKS_PER_ROUND = 4

  # Eight deterministic four-start templates (offsets in units of the current span)
  start_bank = [
    (2, 6, 10, 14),  # A: classic
    (1, 5, 9, 13),   # B: left-shifted
    (3, 7, 11, 15),  # C: right-shifted
    (2, 4, 8, 12),   # D: compressed
    (4, 8, 12, 16),  # E: stretched right
    (2, 7, 10, 15),  # F: offset mix
    (1, 6, 9, 14),   # G: alt-left
    (3, 5, 11, 13),  # H: alt-compact
  ]

  # Seed: a single unit interval (keeps total size within budget for MAIN_ROUNDS = 6)
  T = [(seed_lo, seed_lo + 1.0 * seed_scale)]

  # Helper: build one main-phase round given current T and chosen starts
  def main_round(T, starts, parity_flip):
    lo = min(l for l, r in T)
    hi = max(r for l, r in T)
    delta = hi - lo if hi > lo else 1.0

    # Build the translated blocks, alternating inner order to perturb reuse
    blocks = []
    for b_idx, s in enumerate(starts):
      base = T if ((parity_flip + b_idx) % 2 == 0) else list(reversed(T))
      shift = s * delta - lo
      blocks.append([(l + shift, r + shift) for (l, r) in base])

    # Emit blocks in a round-robin interleaving order every two rounds to vary coupling
    # For even parity_flip: forward interleave, for odd: reverse interleave.
    order = list(range(len(blocks)))
    if parity_flip % 2 == 1:
      order = order[::-1]

    S = []
    maxlen = max(len(b) for b in blocks)
    for i in range(maxlen):
      for idx in order:
        blk = blocks[idx]
        if i < len(blk):
          S.append(blk[i])

    # Four global connectors (Figure-4 style) to propagate colors without creating large cliques
    s0, s1, s2, s3 = starts
    connectors = [
      ((s0 - 1) * delta, (s1 - 1) * delta),  # left cap
      ((s2 + 2) * delta, (s3 + 2) * delta),  # right cap
      ((s0 + 2) * delta, (s2 - 1) * delta),  # cross 1
      ((s1 + 2) * delta, (s3 - 1) * delta),  # cross 2
    ]
    S.extend(connectors)
    return S

  # Main phase: six rounds with eight-pattern cycling and alternating parity every two rounds
  for r in range(MAIN_ROUNDS):
    starts = start_bank[r % len(start_bank)]
    parity_flip = (r // 2)  # alternate interleaving direction every two rounds
    T = main_round(T, starts, parity_flip)

  # Early exit if we've reached the interval budget
  if len(T) >= MAX_TOTAL:
    return T[:MAX_TOTAL]

  # Secondary delta2-driven micro-phase: run a couple of mini-rounds if there is room
  # Each micro-round uses delta2 = span/2 and replicates only a thin sample of T
  # (keeps budget and omega in check while increasing cross-scale FF pressure).
  def thin_sample(U, target_size):
    if not U:
      return []
    n = len(U)
    k = max(8, min(target_size, n))
    step = max(1, n // k)
    sample = U[::step]
    # Trim to exactly k if needed
    if len(sample) > k:
      sample = sample[:k]
    return sample

  # Micro start patterns (deterministic, staggered)
  micro_patterns = [
    (1, 3, 5, 7),
    (2, 4, 6, 8),
    (1, 4, 7, 10),
  ]

  # Leave buffer for final size; plan micro size ≈ up to 30% of remaining budget
  remaining = MAX_TOTAL - len(T)
  per_micro_budget = max(0, int(0.3 * remaining))

  # Determine a conservative sample size per micro block
  # We'll use up to min(512, len(T)//16) intervals per block.
  per_block_sample = 0
  if len(T) > 0:
    per_block_sample = min(512, max(16, len(T) // 16))

  for mr in range(MICRO_ROUNDS):
    if len(T) >= MAX_TOTAL - 64:
      break  # safety
    lo = min(l for l, r in T)
    hi = max(r for l, r in T)
    span = hi - lo if hi > lo else 1.0
    delta2 = max(1.0, span / 2.0)

    starts = micro_patterns[mr % len(micro_patterns)]
    # Build micro blocks from a thin sample of T
    sample = thin_sample(T, per_block_sample)
    micro_blocks = []
    for b_idx in range(MICRO_BLOCKS_PER_ROUND):
      s = starts[b_idx % len(starts)]
      shift = s * delta2 - lo
      base = sample if ((mr + b_idx) % 2 == 0) else list(reversed(sample))
      micro_blocks.append([(l + shift, r + shift) for (l, r) in base])

    # Emit micro blocks in round-robin to maximize mixing
    interleaved = []
    maxlen = max((len(b) for b in micro_blocks), default=0)
    for i in range(maxlen):
      for b in micro_blocks:
        if i < len(b):
          interleaved.append(b[i])

    # Add a small set of micro-connectors (scaled-down versions) to tie the micro blocks
    if len(starts) >= 4:
      s0, s1, s2, s3 = starts[:4]
      micro_connectors = [
        ((s0 + 0.5) * delta2, (s1 + 0.5) * delta2),
        ((s2 - 0.4) * delta2, (s3 - 0.4) * delta2),
        ((s0 + 0.8) * delta2, (s2 - 0.2) * delta2),
        ((s1 + 0.8) * delta2, (s3 - 0.2) * delta2),
      ]
    else:
      micro_connectors = []

    # Prepend micro blocks to act as color blockers for the subsequent main-structure intervals
    T = interleaved + micro_connectors + T

    # Enforce the global cap
    if len(T) > MAX_TOTAL:
      T = T[:MAX_TOTAL]
      break

  # Final tiny sprinkle: a few sparse caps across the full span to lightly bind colors
  if len(T) + 6 <= MAX_TOTAL and T:
    lo = min(l for l, r in T)
    hi = max(r for l, r in T)
    span = max(1.0, hi - lo)
    tiny = []
    for frac in (0.18, 0.37, 0.63, 0.82):
      c = lo + frac * span
      w = 0.04 * span
      tiny.append((c - 0.5 * w, c + 0.5 * w))
    T.extend(tiny)

  # Ensure we never exceed the interval cap
  if len(T) > MAX_TOTAL:
    T = T[:MAX_TOTAL]

  return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()