# EVOLVE-BLOCK-START

def construct_intervals(rounds=3):
  """
  Construct a sequence of intervals presented to FirstFit to push color usage
  while keeping the clique number small, using alternating block order,
  thin bridges, and embedded micro-gadgets with rotated starts.

  Inspired by recursive gadgets in https://arxiv.org/abs/1506.00192.
  """
  # Low-omega spine: multiple disjoint unit intervals
  T = [(0, 1), (2, 3), (4, 5), (6, 7)]

  # Small bank of four-interval gadgets (variations)
  gadget_templates = [
    [(1.0, 5.0), (12.0, 16.0), (4.0, 9.0), (8.0, 13.0)],   # A
    [(0.5, 4.5), (11.0, 15.0), (3.5, 8.5), (7.0, 12.0)],   # B
    [(1.0, 4.0), (6.0, 9.0), (3.0, 7.0), (9.0, 13.0)],     # C
  ]

  # Rotate compact start patterns across rounds (keeps growth modest)
  start_patterns = [
    (2, 6, 10, 14),
    (1, 5, 9, 13),
    (3, 7, 11, 15),
  ]

  for round_idx in range(rounds):
    lo = min(l for l, r in T)
    hi = max(r for l, r in T)
    delta = hi - lo if hi > lo else 1.0
    S = []

    starts = list(start_patterns[round_idx % len(start_patterns)])
    core_gadget = gadget_templates[round_idx % len(gadget_templates)]

    # Alternating inner order and thin bridges (plus occasional second-next bridge)
    for idx, start in enumerate(starts):
      base = T if ((idx + round_idx) % 2 == 0) else list(reversed(T))
      S += [(delta * start + l - lo, delta * start + r - lo) for (l, r) in base]

      if idx + 1 < len(starts):
        nxt = starts[idx + 1]
        S.append((delta * (start + 0.45), delta * (nxt + 0.55)))
      if (round_idx % 2 == 1) and (idx + 2 < len(starts)):
        nxt2 = starts[idx + 2]
        S.append((delta * (start + 0.2), delta * (nxt2 + 0.8)))

    # Core gadget at global scale
    S += [(delta * a, delta * b) for (a, b) in core_gadget]

    # Embedded micro-gadget at 1/16 scale anchored in the middle block
    anchor = starts[len(starts) // 2]
    micro_template = gadget_templates[(round_idx + 1) % len(gadget_templates)]
    micro_scale = 1.0 / 16.0
    S += [(delta * (anchor + a * micro_scale), delta * (anchor + b * micro_scale)) for (a, b) in micro_template]

    # Two sparse cap families to gently increase FF pressure
    max_start = max(starts)
    step1 = max(2, int(max_start // 4))
    step2 = max(3, int(max_start // 3))
    for j in range(2, max_start, step1):
      S.append((delta * (j - 0.6), delta * (j + 0.6)))
    for j in range(3, max_start, step2):
      S.append((delta * (j - 0.2), delta * (j + 1.8)))

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