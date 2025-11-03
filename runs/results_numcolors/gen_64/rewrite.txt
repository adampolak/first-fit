# EVOLVE-BLOCK-START

def construct_intervals(rounds=3):
  """
  Fractal-like multi-scale wave construction.
  Produces a deterministic sequence of intervals (open) for FirstFit.

  Strategy overview:
  - Seed with four disjoint unit intervals to establish a small omega spine.
  - In each round, generate multiple translated copies (blocks) of the current T
    at fixed offsets (offsets grow with round) to create multi-wave overlap.
  - Interleave blocks to promote color diversity; add short bridging intervals
    and a gadget set to pressure color usage without inflating omega.
  - Normalize to non-negative integers and clamp to a reasonable maximum.
  Returns:
    List of (l, r) pairs representing open intervals.
  """
  # Initial spine: four disjoint unit intervals
  T = [(0, 1), (2, 3), (4, 5), (6, 7)]

  for round_idx in range(rounds):
    lo = min(l for l, r in T)
    hi = max(r for l, r in T)
    delta = hi - lo
    if delta <= 0:
      delta = 1

    S = []

    # Offsets for block translations
    if round_idx == 0:
      current_offsets = [2, 6, 10, 14]  # compact first round
    else:
      # Expanded set in later rounds to inject more waves (deterministic)
      current_offsets = [2, 6, 10, 14, 18, 22, 26, 30]

    # Build translated copies (blocks)
    for off in current_offsets:
      translated = [(delta * off + l - lo, delta * off + r - lo) for (l, r) in T]
      S.extend(translated)

    # Gadget intervals to push color usage while keeping omega small
    gadget = [
      (delta * 1, delta * 5),
      (delta * 12, delta * 16),
      (delta * 4, delta * 9),
      (delta * 8, delta * 13),
    ]
    S.extend(gadget)

    # Short bridges between consecutive blocks to couple colors
    for i in range(len(current_offsets) - 1):
      a = delta * (current_offsets[i] + 0.5) - lo
      b = delta * (current_offsets[i + 1] + 0.5) - lo
      if b > a:
        S.append((a, b))

    T = S

  # Normalize to non-negative integers
  if not T:
    return []

  min_l = min(l for l, r in T)
  if min_l < 0:
    T = [(l - min_l, r - min_l) for l, r in T]

  intervals = []
  for (l, r) in T:
    li = int(round(l))
    ri = int(round(r))
    if ri <= li:
      ri = li + 1
    intervals.append((li, ri))

  # Keep within a reasonable bound for the task (less than 10k)
  MAX_INTERVALS = 10000
  if len(intervals) > MAX_INTERVALS:
    intervals = intervals[:MAX_INTERVALS]

  return intervals

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()