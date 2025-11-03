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

  # define a hybrid ordering to diversify interval presentation
  def block_hybrid_order(seq):
    n = len(seq)
    bs = (n + 3) // 4
    blocks = [seq[i*bs:(i+1)*bs] for i in range(4)]
    b0 = blocks[0]
    b1 = list(reversed(blocks[1]))
    b2 = sorted(blocks[2], key=lambda iv: (iv[1]-iv[0], iv[0]))
    b3 = sorted(blocks[3], key=lambda iv: (-(iv[1]-iv[0]), iv[0]))
    return b0 + b1 + b2 + b3

  T = [(0, 1)]
  for _ in range(4):
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

  # --------------------
  # Evaluate multiple arrival orders and pick the one that makes FirstFit use the most colors.
  # This is a targeted, low-risk improvement: same interval set, smarter order selection.
  # --------------------

  # exact FirstFit simulator for open intervals
  def firstfit_colors(seq):
    colors = []  # list of lists of intervals
    for (l, r) in seq:
      placed = False
      for col in colors:
        conflict = False
        for (cl, cr) in col:
          # open intervals overlap iff cl < r and l < cr
          if cl < r and l < cr:
            conflict = True
            break
        if not conflict:
          col.append((l, r))
          placed = True
          break
      if not placed:
        colors.append([(l, r)])
    return len(colors)

  # A small portfolio of orderers (diverse heuristics)
  def order_fractal(seq):
    return list(seq)

  def order_reversed(seq):
    return list(reversed(seq))

  def order_left_first(seq):
    return sorted(seq, key=lambda x: (x[0], x[1]))

  def order_right_first(seq):
    return sorted(seq, key=lambda x: (-x[0], -x[1]))

  def order_short_first(seq):
    return sorted(seq, key=lambda x: (x[1] - x[0], x[0]))

  def order_long_first(seq):
    return sorted(seq, key=lambda x: (-(x[1] - x[0]), x[0]))

  def order_alt_short_long(seq):
    # alternate smallest and largest by length
    by_len = sorted(seq, key=lambda x: (x[1] - x[0], x[0]))
    l, r = 0, len(by_len) - 1
    out = []
    take_left = True
    while l <= r:
      if take_left:
        out.append(by_len[l]); l += 1
      else:
        out.append(by_len[r]); r -= 1
      take_left = not take_left
    return out

  def order_high_degree(seq):
    # place intervals in decreasing order of overlap-degree (greedy heuristic)
    n = len(seq)
    degs = [0] * n
    for i in range(n):
      li, ri = seq[i]
      for j in range(n):
        if i == j:
          continue
        lj, rj = seq[j]
        if li < rj and lj < ri:
          degs[i] += 1
    idxs = sorted(range(n), key=lambda i: (-degs[i], seq[i][1] - seq[i][0]))
    return [seq[i] for i in idxs]

  orderers = [
    ("fractal", order_fractal),
    ("reversed", order_reversed),
    ("left_first", order_left_first),
    ("right_first", order_right_first),
    ("short_first", order_short_first),
    ("long_first", order_long_first),
    ("block_hybrid", block_hybrid_order),
    ("alt_short_long", order_alt_short_long),
    ("high_degree", order_high_degree),
  ]

  best_seq = None
  best_colors = -1
  for name, ord_fn in orderers:
    seq_ord = ord_fn(list(T))
    colors = firstfit_colors(seq_ord)
    # prefer higher FF color count; tie-breaker: smaller sequence (stable) then earlier name
    if colors > best_colors or (colors == best_colors and (best_seq is None or len(seq_ord) < len(best_seq))):
      best_colors = colors
      best_seq = seq_ord

  return best_seq

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