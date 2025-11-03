# EVOLVE-BLOCK-START

def construct_intervals(rounds=2):
  """
  Construct a sequence of intervals presented to FirstFit to push color usage
  while keeping the clique number small, using alternating block order,
  thin bridges, and an embedded micro-gadget.

  Inspired by recursive gadgets in https://arxiv.org/abs/1506.00192.
  """
  T = [(0, 1)]
  for round_idx in range(rounds):
    lo = min(l for l, r in T)
    hi = max(r for l, r in T)
    delta = hi - lo if hi > lo else 1.0
    S = []

    # Compact starts in early round; slightly sparser, offset pattern later to limit blow-up
    if round_idx == 0:
      starts = [2, 4, 6, 8, 10, 12, 14, 16]
    else:
      starts = [2, 4, 6, 8, 10, 12, 15, 18, 22, 26]

    # Alternating inner order and thin bridges between consecutive blocks
    for idx, start in enumerate(starts):
      base = T if (idx % 2 == 0) else list(reversed(T))
      S += [(delta * start + l - lo, delta * start + r - lo) for (l, r) in base]
      if idx + 1 < len(starts):
        nxt = starts[idx + 1]
        # bridge starts/ends well inside blocks to couple colors without spiking omega
        S.append((delta * (start + 0.6), delta * (nxt + 0.4)))

    # Core 4-interval gadget (global scale)
    gadget = [(1.0, 5.0), (12.0, 16.0), (4.0, 9.0), (8.0, 13.0)]
    S += [(delta * a, delta * b) for (a, b) in gadget]

    # Embedded micro-gadget near the median start to frustrate color reuse
    mid = starts[len(starts) // 2]
    micro_scale = 0.25
    S += [(delta * (mid + a * micro_scale), delta * (mid + b * micro_scale)) for (a, b) in gadget]

    # Sparse caps to add gentle pressure without increasing clique size much
    max_start = max(starts)
    step = max(3, max_start // 5)
    for j in range(2, max_start, step):
      S.append((delta * (j - 0.5), delta * (j + 1.2)))

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