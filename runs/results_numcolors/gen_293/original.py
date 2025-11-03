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

  # Six KT-style rounds with deterministic parity interleaving
  for round_idx in range(6):
    starts = template_bank[round_idx % 4]
    lo = min(l for l, r in T)
    hi = max(r for l, r in T)
    delta = hi - lo if hi > lo else 1
    # build translated blocks
    blocks = []
    for start in starts:
      blocks.append([(delta * start + l - lo, delta * start + r - lo) for l, r in T])
    S = []
    # interleave blocks on even rounds, sequential on odd rounds
    if round_idx % 2 == 0:
      maxlen = max(len(b) for b in blocks)
      for i in range(maxlen):
        for blk in blocks:
          if i < len(blk):
            S.append(blk[i])
    else:
      for blk in blocks:
        S.extend(blk)
    # connectors based on the active starts
    s0, s1, s2, s3 = starts
    connectors = [
      (delta * (s0 - 1), delta * (s1 - 1)),  # left cap
      (delta * (s2 + 2), delta * (s3 + 2)),  # right cap
      (delta * (s0 + 2), delta * (s2 - 1)),  # cross 1
      (delta * (s1 + 2), delta * (s3 - 1)),  # cross 2
    ]
    for a, b in connectors:
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

  # Micro-phase B: a thin delta2 micro-round (capacity-guarded) to increase FF colors without raising omega much
  if len(T) < CAP - 16:
    glo = min(l for l, r in T)
    ghi = max(r for l, r in T)
    G = ghi - glo if ghi > glo else 1

    # Thin seed of current T (small and evenly spaced)
    seed_sz = max(8, min(32, len(T) // 300))
    stride = max(1, len(T) // max(1, seed_sz))
    U = [T[i] for i in range(0, len(T), stride)][:seed_sz]

    if U:
      ulo = min(l for l, r in U)

      # Four windows across the global span (stay inside to avoid omega spikes)
      window_fracs = [(0.12, 0.22), (0.35, 0.45), (0.58, 0.68), (0.80, 0.90)]
      micro_blocks = []
      for (fa, fb) in window_fracs:
        win_lo = glo + int(round(fa * G))
        base = win_lo - ulo
        block = [(l + base, r + base) for (l, r) in U]
        micro_blocks.append(block)

      # Interleave the micro-blocks
      micro = []
      maxlen = max(len(b) for b in micro_blocks)
      for i in range(maxlen):
        for blk in micro_blocks:
          if i < len(blk):
            micro.append(blk[i])

      # Fractional-span connectors at the delta2 scale
      micro_connectors = [
        (glo + int(round(0.08 * G)), glo + int(round(0.30 * G))),
        (glo + int(round(0.60 * G)), glo + int(round(0.92 * G))),
        (glo + int(round(0.26 * G)), glo + int(round(0.56 * G))),
        (glo + int(round(0.44 * G)), glo + int(round(0.78 * G))),
      ]
      for a, b in micro_connectors:
        if b > a:
          micro.append((a, b))

      # Capacity guard
      room = CAP - len(T)
      if room > 0:
        if len(micro) > room:
          micro = micro[:room]
        T.extend(micro)

  # Micro-phase C: second guarded micro-round using distinct windows and pin shrinking
  if len(T) < CAP - 8:
    glo2 = min(l for l, r in T)
    ghi2 = max(r for l, r in T)
    G2 = ghi2 - glo2 if ghi2 > glo2 else 1

    # Thin seed again (slightly smaller to limit omega growth)
    seed_sz2 = max(8, min(24, len(T) // 400))
    stride2 = max(1, len(T) // max(1, seed_sz2))
    U2 = [T[i] for i in range(0, len(T), stride2)][:seed_sz2]

    if U2:
      ulo2 = min(l for l, r in U2)

      # New window set to avoid stacking with phase B windows
      windows2 = [(0.05, 0.15), (0.28, 0.38), (0.60, 0.70), (0.82, 0.92)]

      # Short pins to reduce clique growth; slight stagger to avoid co-location
      eps = max(1, G2 // 512)
      blocks2 = []
      for (fa, fb) in windows2:
        win_lo2 = glo2 + int(round(fa * G2))
        base2 = win_lo2 - ulo2
        block = []
        for idx, (l, r) in enumerate(U2):
          mid = (l + r) // 2
          L = mid + base2 - (eps // 2) + (idx % 3)
          R = L + eps
          if R > L:
            block.append((L, R))
        blocks2.append(block)

      # Interleave with a different order to diversify coupling
      micro2 = []
      if blocks2:
        maxlen2 = max(len(b) for b in blocks2)
        order2 = [3, 1, 0, 2]
        for i in range(maxlen2):
          for idx in order2:
            blk = blocks2[idx]
            if i < len(blk):
              micro2.append(blk[i])

      # Very short cap-pins to tie windows without inflating omega
      caps2 = [
        (glo2 + int(round(0.09 * G2)), glo2 + int(round(0.09 * G2)) + eps),
        (glo2 + int(round(0.64 * G2)), glo2 + int(round(0.64 * G2)) + eps),
      ]
      for a, b in caps2:
        if b > a:
          micro2.append((a, b))

      # Capacity guard
      room2 = CAP - len(T)
      if room2 > 0:
        if len(micro2) > room2:
          micro2 = micro2[:room2]
        T.extend(micro2)

  # Micro-phase D: second four-window micro-phase with distinct windows and long-range connectors
  if len(T) < CAP - 4:
    glo3 = min(l for l, r in T)
    ghi3 = max(r for l, r in T)
    G3 = ghi3 - glo3 if ghi3 > glo3 else 1

    # Thin seed (moderate) from current T for stronger mixing
    seed_sz3 = max(8, min(36, len(T) // 300))
    stride3 = max(1, len(T) // max(1, seed_sz3))
    U3 = [T[i] for i in range(0, len(T), stride3)][:seed_sz3]

    if U3:
      ulo3 = min(l for l, r in U3)

      # Distinct window pattern to complement earlier phases
      window_fracs3 = [(0.09, 0.19), (0.32, 0.42), (0.64, 0.74), (0.86, 0.94)]
      blocks3 = []
      for (fa, fb) in window_fracs3:
        win_lo3 = glo3 + int(round(fa * G3))
        base3 = win_lo3 - ulo3
        block = [(l + base3, r + base3) for (l, r) in U3]
        # Reverse every second block to vary internal order
        if int(round(fa * 100)) % 2 == 0:
          block = list(reversed(block))
        blocks3.append(block)

      # Interleave using a nontrivial order to break symmetry
      micro3 = []
      if blocks3:
        maxlen3 = max(len(b) for b in blocks3)
        order3 = [0, 2, 1, 3]
        for i in range(maxlen3):
          for idx in order3:
            blk = blocks3[idx]
            if i < len(blk):
              micro3.append(blk[i])

      # Long-range, pairwise-disjoint connectors to avoid omega spikes
      connectors3 = [
        (glo3 + int(round(0.06 * G3)), glo3 + int(round(0.28 * G3))),
        (glo3 + int(round(0.41 * G3)), glo3 + int(round(0.59 * G3))),
        (glo3 + int(round(0.72 * G3)), glo3 + int(round(0.95 * G3))),
      ]
      for a, b in connectors3:
        if b > a:
          micro3.append((a, b))

      # Capacity guard
      room3 = CAP - len(T)
      if room3 > 0:
        if len(micro3) > room3:
          micro3 = micro3[:room3]
        T.extend(micro3)

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