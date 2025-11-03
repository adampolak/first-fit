# EVOLVE-BLOCK-START

def construct_intervals(rounds=4):
  """
  Construct a sequence of intervals in presentation order with diversified
  start-patterns, connector templates, light bridges and caps to increase
  FirstFit color pressure while keeping the clique number small.

  Changes vs. the previous simple replicator:
  - Seed with four disjoint unit intervals to create richer early interactions.
  - Cycle among several deterministic start-patterns to disrupt regular reuse.
  - Cycle a small bank of 4-interval connector templates (one per round).
  - Add light bridges and short caps between translated blocks to couple colors.
  - Default rounds increased to 4; rounds is configurable.
  """
  # Seed: several disjoint unit intervals to keep omega low but enable coupling
  T = [
    (0.0, 1.0),
    (2.0, 3.0),
    (4.0, 5.0),
    (6.0, 7.0)
  ]

  # Cycle of start-patterns (rotating choices to avoid predictable reuse)
  start_patterns = [
      (2, 6, 10, 14),
      (1, 5, 9, 13),
      (3, 7, 11, 15),
      (2, 4, 8, 12),
  ]

  # Extended rotation sequence to diversify start pattern usage across rounds
  rotation_seq = [0, 1, 2, 3, 1, 2, 3, 0]

  # Small bank of 4-interval connector templates (keeps per-round additive cost small)
  template_bank = [
      ((1, 5), (12, 16), (4, 9), (8, 13)),
      ((0.5, 4.5), (11, 15), (3.5, 8.5), (7, 12)),
      ((1, 4), (6, 9), (3, 7), (9, 13)),
      ((2, 6), (7, 11), (0, 3), (10, 14)),
  ]

  for round_idx in range(rounds):
    lo = min(l for l, r in T)
    hi = max(r for l, r in T)
    span = hi - lo if hi > lo else 1.0

    # pick start pattern and template using extended rotation sequence
    rot_idx = rotation_seq[round_idx % len(rotation_seq)]
    starts = list(start_patterns[rot_idx])
    template = template_bank[rot_idx]

    S = []

    # vary block insertion order each round to perturb arrival ordering
    if round_idx % 3 == 0:
      block_order = starts
    elif round_idx % 3 == 1:
      block_order = list(reversed(starts))
    else:
      block_order = starts[1:] + starts[:1]

    # Build translated copies with alternating inner order for T
    for i, st in enumerate(block_order):
      base = T if ((round_idx + i) % 2 == 0) else list(reversed(T))
      base_off = span * st - lo
      # append clones (preserving presentation order)
      for (l, r) in base:
        S.append((l + base_off, r + base_off))

      # light bridge: a short interval that overlaps the tail of this block and head of the next
      if i + 1 < len(block_order):
        nxt = block_order[i + 1]
        bstart = span * (st + 0.5) - lo
        bend = span * (nxt + 0.5) - lo
        # make it short compared to span to avoid increasing omega much
        S.append((min(bstart, bend) + 0.1 * span * 1e-6, max(bstart, bend) - 0.1 * span * 1e-6))

    # add the 4-interval connector gadget (scaled by span)
    for (a, b) in template:
      S.append((span * a, span * b))

    # add a few short caps sparsely between blocks to press colors without raising omega
    max_st = max(starts)
    step = max(2, int(max_st // 3))
    for j in range(2, max_st, step):
      a = span * (j - 0.4)
      b = span * (j + 0.6)
      # make these caps short relative to span
      S.append((a + 0.01 * span, b - 0.01 * span))

    # add long spines across full block range every other round to couple distant colors
    if round_idx % 2 == 1:
      left = min(starts) - 0.3
      right = max(starts) + 0.3
      S.append((span * left, span * right))
      S.append((span * (left + 0.5), span * (right + 0.5)))
    # update T for next round
    T = S

  return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()