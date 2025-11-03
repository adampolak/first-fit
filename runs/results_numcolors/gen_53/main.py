# EVOLVE-BLOCK-START

def construct_intervals(rounds=6, rotate_starts=True, reverse_block_parity=True, interleave_blocks=True, phase2_iters=1):
  """
  Deterministic rotated four-block Kierstead–Trotter style expansion with block interleaving.

  Parameters:
    rounds (int): main expansion depth; 6 yields ~9556 intervals (near cap).
    rotate_starts (bool): rotate the four translated starts across rounds to disrupt
                          repeating overlap patterns.
    reverse_block_parity (bool): reverse the order of T for every odd block within a round.
    interleave_blocks (bool): interleave translated blocks round-robin to enhance color mixing.
    phase2_iters (int): optional micro-scale follow-up iterations (kept 1 by default).

  Returns:
    intervals: list of (l, r) integer tuples, open intervals, in FF presentation order.
  """

  # Rotation cycle for the four translated copies per round.
  # These are carefully chosen to preserve coupling while varying interactions:
  start_patterns = [
    [2, 6, 10, 14],  # classic
    [1, 5, 9, 13],   # left-shifted
    [3, 7, 11, 15],  # right-shifted
    [2, 4, 8, 12],   # compressed left pair
  ]

  # Base seed (unit interval). Using 1 seed keeps growth within limit at 6 rounds.
  T = [(0, 1)]

  # Main deterministic rotated expansion
  rounds = max(1, int(rounds))
  for round_idx in range(rounds):
    # Span of the current set
    lo = min(l for l, r in T)
    hi = max(r for l, r in T)
    delta = hi - lo  # integer span

    # Choose the starts for this round
    if rotate_starts:
      starts = start_patterns[round_idx % len(start_patterns)]
    else:
      starts = start_patterns[0]

    # Build four translated blocks with optional reversal for parity
    blocks = []
    for b_idx, s in enumerate(starts):
      block_src = T[::-1] if (reverse_block_parity and (b_idx % 2 == 1)) else T
      base = s * delta - lo
      blocks.append([(l + base, r + base) for (l, r) in block_src])

    # Interleave blocks round-robin to enhance color mixing if enabled
    if interleave_blocks:
      S = []
      maxlen = max(len(b) for b in blocks)
      order = list(range(len(blocks)))
      # rotate block order by round index to disrupt patterns
      rot = round_idx % len(order)
      order = order[rot:] + order[:rot]
      for i in range(maxlen):
        for idx in order:
          if i < len(blocks[idx]):
            S.append(blocks[idx][i])
    else:
      S = []
      for blk in blocks:
        S.extend(blk)

    # Deterministically computed connectors derived from the selected starts.
    # Generic formulas that recover Figure 4 at starts=[2,6,10,14]:
    s0, s1, s2, s3 = starts
    connectors = [
      ( (s0 - 1) * delta, (s1 - 1) * delta ),  # left cap
      ( (s2 + 2) * delta, (s3 + 2) * delta ),  # right cap
      ( (s0 + 2) * delta, (s2 - 1) * delta ),  # cross 1
      ( (s1 + 2) * delta, (s3 - 1) * delta ),  # cross 2
    ]
    # Append connectors in a color-chaining order
    for (a, b) in connectors:
      S.append((a, b))

    T = S

  # Optional tiny second phase (disabled by default to keep count near 9556)
  # This can add sparse, long-range gadgets at a smaller scale without
  # blowing up n or omega. Left here as deterministic, but phase2_iters=0 by default.
  if phase2_iters > 0:
    for k in range(phase2_iters):
      lo = min(l for l, r in T)
      hi = max(r for l, r in T)
      delta = hi - lo
      # Use a smaller subscale to sprinkle a few caps; keep them sparse.
      d2 = max(1, delta // 4)
      # Choose a rotating micro-template; keep count constant per phase.
      micro = [
        (lo + 1 * d2, lo + 5 * d2),
        (hi - 6 * d2, hi - 2 * d2),
        (lo + 3 * d2, lo + 8 * d2),
        (hi - 8 * d2, hi - 3 * d2),
      ]
      # Order interleaved with existing to increase FF pressure slightly.
      T.extend(micro)

  return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()