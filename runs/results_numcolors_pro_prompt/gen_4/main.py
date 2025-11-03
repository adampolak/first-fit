# EVOLVE-BLOCK-START

def construct_intervals():
  """
  Omega-guarded, chained multi-block expansion to improve FirstFit ratio
  while keeping offline optimum (omega) bounded by 10.
  We generalize the 4-block gadget to m blocks per round with sliding connectors
  that couple colors across all blocks. We search over several widths (m) and
  up to 4 rounds and choose the largest valid construction (proxy for higher FF),
  with a strict omega cap and size bound. A robust fallback preserves the
  original two-round construction if needed.
  Returns:
    intervals: list of tuples (l, r) open intervals
  """

  def compute_omega(intervals):
    # compute maximum number of intervals covering any point
    events = []
    for a, b in intervals:
      events.append((a, 1))   # start
      events.append((b, -1))  # end
    # For open intervals, end before start at the same coordinate
    events.sort(key=lambda x: (x[0], x[1]))
    current = 0
    best = 0
    for _, typ in events:
      if typ == -1:
        current -= 1
      else:
        current += 1
        if current > best:
          best = current
    return best

  def expand_once(T, m):
    # Replicate T at m shifted positions and add sliding connector families
    lo = min(l for l, r in T)
    hi = max(r for l, r in T)
    delta = hi - lo
    S = []
    # m replicated blocks at starts 2,6,10,..., (4m-2)
    for i in range(m):
      start = 4 * i + 2
      S += [(delta * start + l - lo, delta * start + r - lo) for l, r in T]
    # Sliding connectors: four short families that couple adjacent blocks
    # Family A: (4i+1, 4i+5), i = 0..m-2
    for i in range(max(0, m - 1)):
      S.append((delta * (4 * i + 1), delta * (4 * i + 5)))
    # Family B: (4i+4, 4i+9), i = 0..m-3
    for i in range(max(0, m - 2)):
      S.append((delta * (4 * i + 4), delta * (4 * i + 9)))
    # Family C and D only when m >= 4 to avoid overshooting omega early
    if m >= 4:
      # Family C: (4i+8, 4i+13), i = 0..m-4
      for i in range(max(0, m - 3)):
        S.append((delta * (4 * i + 8), delta * (4 * i + 13)))
      # Family D: (4i+12, 4i+16), i = 0..m-4
      for i in range(max(0, m - 3)):
        S.append((delta * (4 * i + 12), delta * (4 * i + 16)))
    return S

  # Search over widths and rounds under omega<=10 and size constraint
  candidate_widths = [9, 8, 7, 6]
  max_rounds = 4
  max_n = 6000  # keep constructions compact for better combined score
  best_T = None
  best_w = None

  for m in candidate_widths:
    T = [(0, 1)]
    local_best_T = None
    local_best_w = None
    for _ in range(max_rounds):
      T = expand_once(T, m)
      if len(T) > max_n:
        break
      w = compute_omega(T)
      if w <= 10:
        local_best_T = list(T)
        local_best_w = w
      else:
        break
    # prefer larger sequences under the omega cap; tie-break by larger omega (harder)
    if local_best_T is not None:
      if (best_T is None or
          len(local_best_T) > len(best_T) or
          (len(local_best_T) == len(best_T) and (best_w is None or local_best_w > best_w))):
        best_T = local_best_T
        best_w = local_best_w

  if best_T is None:
    # Fall back to the original two-round construction to guarantee omega <= 10
    T = [(0, 1)]
    for _ in range(2):
      lo = min(l for l, r in T)
      hi = max(r for l, r in T)
      delta = hi - lo
      S = []
      for start in (2, 6, 10, 14):
        S += [(delta * start + l - lo, delta * start + r - lo) for l, r in T]
      S += [
        (delta * 1, delta * 5),
        (delta * 12, delta * 16),
        (delta * 4, delta * 9),
        (delta * 8, delta * 13)
      ]
      T = S
    best_T = T

  return best_T

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