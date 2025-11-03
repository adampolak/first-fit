# EVOLVE-BLOCK-START

def construct_intervals():
  """
  Construct a sequence of intervals of real line,
  in the order in which they are presented to FirstFit,
  so that it maximizes the number of colors used by FirstFit
  divided by the maximum number of intervals that cover a single point

  The initial implementation uses the construction from
  Figure 4 in https://arxiv.org/abs/1506.00192

  Returns:
    intervals: list of tuples, each tuple (l, r) represents an open interval from l to r
  """

  # Cycle through multiple start-pattern templates to break periodic reuse
  start_patterns = [(2, 6, 10, 14), (1, 5, 9, 13), (3, 7, 11, 15), (2, 4, 8, 12)]
  T = [(0, 1)]
  for round_idx in range(6):
    lo = min(l for l, r in T)
    hi = max(r for l, r in T)
    delta = hi - lo
    S = []
    # choose the start tuple for this round
    starts = start_patterns[round_idx % len(start_patterns)]

    # build translated copies with alternating parity to disrupt reuse
    for idx, start in enumerate(starts):
      base_seq = T if (idx % 2 == 0) else list(reversed(T))
      S += [(delta * start + l - lo, delta * start + r - lo) for l, r in base_seq]
      # add a short bridge towards the next block (if any)
      if idx + 1 < len(starts):
        nxt = starts[idx + 1]
        b_lo = delta * (start + 0.35)
        b_hi = delta * (nxt - 0.35)
        if b_hi > b_lo:
          # short interior bridge that overlaps many active colors but avoids large cliques
          S.append((b_lo, min(b_lo + max(1.0, 0.12 * delta), b_hi)))

    # dynamic connectors based on the current starts
    s0, s1, s2, s3 = starts
    connectors = [
      (delta * (s0 - 1), delta * (s1 - 1)),  # left cap
      (delta * (s2 + 2), delta * (s3 + 2)),  # right cap
      (delta * (s0 + 2), delta * (s2 - 1)),  # cross 1
      (delta * (s1 + 2), delta * (s3 - 1)),  # cross 2
    ]
    for a, b in connectors:
      S.append((a, b))

    # add a single sparse spine across the round to couple distant colors
    spine_lo = delta * (s0 - 0.5)
    spine_hi = delta * (s3 + 0.5)
    if spine_hi > spine_lo:
      S.append((spine_lo, spine_hi))

    # tiny caps near each start to add local pressure with minimal omega impact
    for st in starts:
      cap_lo = delta * (st + 0.10)
      cap_hi = cap_lo + max(1.0, 0.10 * delta)
      S.append((cap_lo, cap_hi))

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