# EVOLVE-BLOCK-START

def construct_intervals():
  """
  Construct a sequence of open intervals in the order presented to FirstFit,
  aiming to maximize the FirstFit-to-omega ratio while keeping omega ≤ 10.
  Strategy:
  - Use a low-omega "spine" of long intervals (seeded by disjoint units),
    and recursively build translated waves across 3 rounds.
  - Alternate inner order across translated copies to hinder FirstFit reuse.
  - Add light bridges and sparse caps to couple colors across blocks
    without significantly increasing the clique number.
  - Inject the 4-interval gadget each round (Figure 4-style) to raise pressure.
  Returns:
    list[(l, r)]: open intervals
  """
  # Spine seed: 4 disjoint unit intervals
  T = [(0.0, 1.0), (2.0, 3.0), (4.0, 5.0), (6.0, 7.0)]

  # Start-patterns cycle: compact, then denser, then shifted-compact
  start_patterns = [
    (2, 6, 10, 14),
    (2, 4, 6, 8, 10),
    (3, 7, 11, 15),
  ]

  # Four-interval gadget per round (scaled by delta)
  gadget = [(1, 5), (12, 16), (4, 9), (8, 13)]

  rounds = 3
  for round_idx in range(rounds):
    lo = min(l for l, r in T)
    hi = max(r for l, r in T)
    delta = hi - lo if hi > lo else 1.0

    S = []
    starts = list(start_patterns[round_idx % len(start_patterns)])

    # Build translated copies with alternating inner order
    for idx, start in enumerate(starts):
      base = T if (idx % 2 == 0) else list(reversed(T))
      shift = delta * start
      S.extend((shift + l - lo, shift + r - lo) for (l, r) in base)

      # Light bridge to next start to couple colors; keeps omega modest
      if idx + 1 < len(starts):
        nxt = starts[idx + 1]
        S.append((delta * (start + 0.35), delta * (nxt + 0.65)))

    # Inject local gadget (scaled)
    S.extend((delta * a, delta * b) for (a, b) in gadget)

    # Sparse caps: short spans around interiors of blocks
    step = max(2, int(max(starts) // 4))
    for j in range(2, max(starts), step):
      S.append((delta * (j - 0.5), delta * (j + 0.5)))
    for j in range(3, max(starts), step + 1):
      S.append((delta * (j - 0.2), delta * (j + 1.2)))

    T = S

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