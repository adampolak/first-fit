# EVOLVE-BLOCK-START

def construct_intervals(seed_count=1):
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

  # Capacity guard to keep the total count under 10000 and allow a micro-phase
  CAP = 9800

  # Rotate among four strong start-pattern templates to couple colors across scales
  template_bank = [
    (2, 6, 10, 14),
    (1, 5, 9, 13),
    (3, 7, 11, 15),
    (4, 8, 12, 16),
  ]

  # Seed with multiple disjoint unit intervals if requested (new capability)
  if seed_count <= 1:
    T = [(0, 1)]
  else:
    step = 3
    T = [(i * step, i * step + 1) for i in range(seed_count)]

  # Six KT-style rounds with enhanced parity interleaving and sparse late cross4 connectors
  for round_idx in range(6):
    starts = template_bank[round_idx % 4]
    lo = min(l for l, r in T)
    hi = max(r for l, r in T)
    delta = hi - lo if hi > lo else 1
    # build translated blocks with parity-based internal reversals to disrupt FF reuse
    blocks = []
    for idx, start in enumerate(starts):
      base = delta * start - lo
      src = T[::-1] if ((round_idx + idx) % 2 == 1) else T
      blocks.append([(base + l, base + r) for l, r in src])
    S = []
    # Interleave on all rounds; odd rounds use reversed block order to increase mixing
    maxlen = max(len(b) for b in blocks)
    order = list(range(len(blocks)))
    if round_idx % 2 == 1:
      order.reverse()
    for i in range(maxlen):
      for j in order:
        blk = blocks[j]
        if i < len(blk):
          S.append(blk[i])
    # connectors based on the active starts
    s0, s1, s2, s3 = starts
    connectors = [
      (delta * (s0 - 1), delta * (s1 - 1)),  # left cap
      (delta * (s2 + 2), delta * (s3 + 2)),  # right cap
      (delta * (s0 + 2), delta * (s2 - 1)),  # cross 1
      (delta * (s1 + 2), delta * (s3 - 1)),  # cross 2
    ]
    # Sparse late cross4 connectors to couple distant towers (last two rounds only)
    if round_idx >= 4:
      connectors.extend([
        (delta * (s0 + 1), delta * (s3 - 1)),
        (delta * (s1 - 1), delta * (s2 + 1)),
      ])
    for a, b in connectors:
      if b > a:
        S.append((a, b))
    T = S

  # Micro-phase A: add three long-range caps to boost late FirstFit color usage while controlling omega
  lo = min(l for l, r in T)
  hi = max(r for l, r in T)
  delta = hi - lo if hi > lo else 1
  def cap_at(a_frac, b_frac):
    l0 = lo + max(1, int(round(a_frac * delta)))
    r0 = lo + max(1, int(round(b_frac * delta)))
    return (l0, r0)
  capA = cap_at(0.08, 0.60)
  capB = cap_at(0.25, 0.75)
  capC = cap_at(0.75, 0.92)
  for cap in (capA, capB, capC):
    if cap[1] > cap[0]:
      T.append(cap)

  # Micro-phase B and C: two thin delta2 rounds with distinct windows and interleaving, capacity-guarded
  if len(T) < CAP - 16:
    glo = min(l for l, r in T)
    ghi = max(r for l, r in T)
    G = ghi - glo if ghi > glo else 1

    def build_micro(U, window_fracs, reverse_blocks=False, reverse_order=False):
      if not U:
        return []
      ulo = min(l for l, r in U)
      # Build blocks anchored at windows
      blocks = []
      for (fa, fb) in window_fracs:
        win_lo = glo + int(round(fa * G))
        base = win_lo - ulo
        blk = [(l + base, r + base) for (l, r) in U]
        if reverse_blocks:
          blk = list(reversed(blk))
        blocks.append(blk)
      # Interleave blocks
      micro = []
      maxlen = max(len(b) for b in blocks)
      order = list(range(len(blocks)))
      if reverse_order:
        order.reverse()
      for i in range(maxlen):
        for idx in order:
          b = blocks[idx]
          if i < len(b):
            micro.append(b[i])
      # Fractional-span connectors at global scale
      connectors = [
        (glo + int(round(0.08 * G)), glo + int(round(0.30 * G))),
        (glo + int(round(0.60 * G)), glo + int(round(0.92 * G))),
        (glo + int(round(0.26 * G)), glo + int(round(0.56 * G))),
        (glo + int(round(0.44 * G)), glo + int(round(0.78 * G))),
      ]
      for a, b in connectors:
        if b > a:
          micro.append((a, b))
      return micro

    # Slightly denser thin seed to increase active-color sampling
    seed_sz = max(12, min(40, len(T) // 280))
    stride = max(1, len(T) // max(1, seed_sz))
    U = [T[i] for i in range(0, len(T), stride)][:seed_sz]

    # Micro-phase B: original windows, forward order
    micro_B = []
    if U:
      window_fracs_B = [(0.12, 0.22), (0.35, 0.45), (0.58, 0.68), (0.80, 0.90)]
      micro_B = build_micro(U, window_fracs_B, reverse_blocks=False, reverse_order=False)

    # Apply micro B under capacity
    room = CAP - len(T)
    if room > 0 and micro_B:
      if len(micro_B) > room:
        micro_B = micro_B[:room]
      T.extend(micro_B)

    # Micro-phase C: shifted windows and reversed interleaving to diversify
    if len(T) < CAP - 16:
      # Recompute room and optionally refresh U thinly (reuse stride to keep coherence)
      room = CAP - len(T)
      U2 = [T[i] for i in range(0, len(T), stride)][:seed_sz] if room > 0 else []
      micro_C = []
      if U2:
        window_fracs_C = [(0.18, 0.26), (0.49, 0.57), (0.70, 0.78)]
        micro_C = build_micro(U2, window_fracs_C, reverse_blocks=True, reverse_order=True)
      if micro_C:
        room = CAP - len(T)
        if room > 0:
          if len(micro_C) > room:
            micro_C = micro_C[:room]
          T.extend(micro_C)

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