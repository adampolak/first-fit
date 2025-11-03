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

  # Capacity guard to keep total intervals < 10000
  CAP = 9800

  # Stage 1: Deterministic KT-style spine (no interleaving), 6 rounds, classic connectors.
  spine_starts = (2, 6, 10, 14)

  # Seed: keep a single unit interval; multi-seed tends to inflate omega too early.
  T = [(0, 1)]

  # Perform up to six spine rounds unless capacity would be exceeded (it won't for this growth).
  for ridx in range(6):
    # Compute span and delta
    lo = min(l for l, r in T)
    hi = max(r for l, r in T)
    delta = hi - lo
    if delta <= 0:
      delta = 1

    # Build four translated blocks
    blocks = []
    for s in spine_starts:
      base = s * delta - lo
      block = [(l + base, r + base) for (l, r) in T]
      blocks.append(block)

    # Interleave blocks on even rounds; keep sequential on odd rounds
    S = []
    if ridx % 2 == 0:
      maxlen = max(len(b) for b in blocks)
      for i in range(maxlen):
        for blk in blocks:
          if i < len(blk):
            S.append(blk[i])
    else:
      for blk in blocks:
        S.extend(blk)

    # Classic connectors (Figure 4 style)
    s0, s1, s2, s3 = spine_starts
    connectors = [
      ((s0 - 1) * delta, (s1 - 1) * delta),  # left cap
      ((s2 + 2) * delta, (s3 + 2) * delta),  # right cap
      ((s0 + 2) * delta, (s2 - 1) * delta),  # cross 1
      ((s1 + 2) * delta, (s3 - 1) * delta),  # cross 2
    ]
    S.extend(connectors)
    T = S

    if len(T) > CAP:
      break

  # Return the spine-only construction to avoid omega inflation from micro-phases.
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