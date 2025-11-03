# EVOLVE-BLOCK-START

def construct_intervals(rounds=3,
                        seed_lo=0.0,
                        seed_scale=1.0,
                        rotate_starts=True,
                        reverse_block_parity=True,
                        interleave_blocks=True,
                        phase2_micro=True):
  """
  Improved recursive wave construction.

  Key changes vs. the previous simple replicator:
  - Seed with four disjoint unit intervals to create a richer initial color spine.
  - Rotate among several deterministic 4-start patterns to diversify geometry between rounds.
  - Cycle a small bank of 4-interval connector templates (Figure-4 style).
  - Alternate inner block order (reverse every other block) and vary block insertion order by round.
  - Insert short bridges between translated blocks and sparse short caps.
  - Optional small micro-phase appended at the end to sprinkle intervals that touch many existing blocks.

  Parameters:
    rounds: number of recursive expansion rounds (default 3)
    seed_lo / seed_scale: knobs to shift/scale seed placement
    rotate_starts, reverse_block_parity, interleave_blocks, phase2_micro: toggles for behavior

  Returns:
    intervals: list of tuples (l, r) as open intervals in presentation order.
  """
  # Seed with several disjoint unit intervals (a small "spine" to raise early-color diversity)
  T = [
    (seed_lo + 0.0 * seed_scale, seed_lo + 1.0 * seed_scale),
    (seed_lo + 2.0 * seed_scale, seed_lo + 3.0 * seed_scale),
    (seed_lo + 4.0 * seed_scale, seed_lo + 5.0 * seed_scale),
    (seed_lo + 6.0 * seed_scale, seed_lo + 7.0 * seed_scale)
  ]

  # Deterministic rotation of 4-start patterns (keeps per-round branching small)
  start_patterns = [
    (2, 6, 10, 14),
    (1, 5, 9, 13),
    (3, 7, 11, 15),
    (2, 4, 8, 12),
  ]

  # Small bank of 4-interval connector templates (Figure-4 style gadgets)
  template_bank = [
    ((1, 5), (12, 16), (4, 9), (8, 13)),
    ((0.5, 4.5), (11, 15), (3.5, 8.5), (7, 12)),
    ((1, 4), (6, 9), (3, 7), (9, 13)),
    ((2, 6), (7, 11), (0, 3), (10, 14))
  ]

  for round_idx in range(rounds):
    lo = min(l for l, r in T)
    hi = max(r for l, r in T)
    span = hi - lo if hi > lo else 1.0

    # choose deterministic start-pattern/template for this round
    starts = list(start_patterns[round_idx % len(start_patterns)]) if rotate_starts else list(start_patterns[0])
    template = template_bank[round_idx % len(template_bank)]

    S = []

    # rotate block insertion order to perturb arrival ordering each round
    if round_idx % 4 == 0:
      block_order = starts
    elif round_idx % 4 == 1:
      block_order = list(reversed(starts))
    elif round_idx % 4 == 2:
      block_order = starts[1:] + starts[:1]
    else:
      # another deterministic shuffle
      block_order = starts[2:] + starts[:2]

    # Build translated copies with alternating inner order for T to reduce color reuse
    for i, st in enumerate(block_order):
      base = T if not reverse_block_parity or ((round_idx + i) % 2 == 0) else list(reversed(T))
      shift = span * st - lo
      for (l, r) in base:
        S.append((l + shift, r + shift))

      # short bridge overlapping tail of this block and head of next block to couple colors
      if i + 1 < len(block_order):
        nxt = block_order[i + 1]
        bstart = span * (st + 0.5) - lo
        bend = span * (nxt + 0.5) - lo
        left = min(bstart, bend) + 0.02 * span
        right = max(bstart, bend) - 0.02 * span
        if right > left:
          # keep bridge short relative to span so we don't blow up omega
          S.append((left, right))

    # Add the 4-interval connector gadget (scaled by current span)
    for (a, b) in template:
      S.append((span * a, span * b))

    # Add sparse short caps between blocks to press color usage without creating large cliques
    max_st = max(starts)
    step = max(2, int(max(1, max_st // 3)))
    for j in range(2, max_st, step):
      a = span * (j - 0.35)
      b = span * (j + 0.55)
      # trim slightly to keep them short
      S.append((a + 0.01 * span, b - 0.01 * span))

    # adopt interleaving strategy (round-robin consumption) if requested
    if interleave_blocks:
      # perform a light interleaving pass to further mix the order
      interleaved = []
      block_len = max(1, len(T))
      # slice S into blocks of block_len then emit round-robin
      slices = [S[k:k + block_len] for k in range(0, len(S), block_len)]
      maxslice = max(len(sl) for sl in slices) if slices else 0
      for pos in range(maxslice):
        for sl in slices:
          if pos < len(sl):
            interleaved.append(sl[pos])
      T = interleaved
    else:
      T = S

  # Optional small micro-phase appended at the end:
  # sprinkle a few short intervals that touch many of the existing translated blocks.
  if phase2_micro:
    lo = min(l for l, r in T)
    hi = max(r for l, r in T)
    span = max(1.0, hi - lo)
    micros = []
    # place a few tiny caps at fractional positions that often intersect many blocks
    for frac in (0.22, 0.5, 0.78):
      center = lo + frac * span
      w = 0.06 * span
      micros.append((center - 0.5 * w, center + 0.5 * w))
    # append micro-phase intervals near the end (so they are presented late)
    T.extend(micros)

  return T

  # return [  # Figure 3, OPT=2, FF=4
  #   (2,3),
  #   (6,7),
  #   (10,11),
  #   (14,15),
  #   (1,5),
  #   (12,16),
  #   (4,9),
  #   (8,13),
  # ]

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()