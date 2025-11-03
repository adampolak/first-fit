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

  # Small template bank (rotate per KT round for better cross-scale coupling)
  template_bank = [
    (2, 6, 10, 14),
    (1, 5, 9, 13),
    (3, 7, 11, 15),
    (4, 8, 12, 16),
  ]

  # Start with a single seed interval
  T = [(0, 1)]

  # Improved KT spine: rotate templates and introduce parity interleaving to mix colors.
  for round_idx in range(6):
    lo = min(l for l, r in T)
    hi = max(r for l, r in T)
    delta = hi - lo
    if delta <= 0:
      delta = 1

    # rotate template to diversify coupling
    starts = template_bank[round_idx % len(template_bank)]

    # build translated blocks
    blocks = []
    for start in starts:
      base = start * delta - lo
      blocks.append([(l + base, r + base) for l, r in T])

    S = []
    # parity-based policy: interleave on even rounds, reverse-block sequence on odd rounds
    if round_idx % 2 == 0:
      maxlen = max(len(b) for b in blocks)
      for i in range(maxlen):
        for blk in blocks:
          if i < len(blk):
            S.append(blk[i])
    else:
      for blk in reversed(blocks):
        S.extend(blk)

    # connectors based on the per-round starts (classic KT caps)
    s0, s1, s2, s3 = starts
    connectors = [
      ((s0 - 1) * delta, (s1 - 1) * delta),  # left cap
      ((s2 + 2) * delta, (s3 + 2) * delta),  # right cap
      ((s0 + 2) * delta, (s2 - 1) * delta),  # cross 1
      ((s1 + 2) * delta, (s3 - 1) * delta),  # cross 2
    ]
    for a, b in connectors:
      S.append((a, b))
    T = S

  # Stage 2: dual-window micro-phases (primary and alternate) using thin sampling.
  # This mixes colors across distant regions without creating large cliques.
  if len(T) < CAP - 16:
    def thin_seed(current_T, max_seed):
      n = len(current_T)
      if n == 0 or max_seed <= 0:
        return []
      step = max(1, n // max_seed)
      return current_T[::step][:max_seed]

    def micro_windows_round(current_T, round_id, budget, alt=False):
      if budget <= 0 or not current_T:
        return []
      glo = min(l for l, r in current_T)
      ghi = max(r for l, r in current_T)
      G = max(1, ghi - glo)

      # thin deterministic seed
      seed_sz = max(8, min(40, len(current_T) // (220 if alt else 300)))
      stride = max(1, len(current_T) // max(1, seed_sz))
      U = [current_T[i] for i in range(0, len(current_T), stride)][:seed_sz]
      if not U:
        return []

      ulo = min(l for l, r in U)

      # window families: primary has small deterministic shift, alternate uses slightly offset windows
      if not alt:
        shift = (round_id % 3) * 0.02
        window_fracs = [
          (0.12 + shift, 0.22 + shift),
          (0.35 + shift, 0.45 + shift),
          (0.58 + shift, 0.68 + shift),
          (0.80 + shift, 0.90 + shift),
        ]
      else:
        window_fracs = [
          (0.05, 0.15),
          (0.28, 0.38),
          (0.60, 0.70),
          (0.82, 0.92),
        ]

      # clamp windows
      window_fracs = [(max(0.03, a), min(0.97, b)) for (a, b) in window_fracs]

      # build blocks aligned to windows
      blocks = []
      for idx, (fa, fb) in enumerate(window_fracs):
        win_lo = glo + int(round(fa * G))
        base = win_lo - ulo
        block = [(l + base, r + base) for (l, r) in U]
        # alternate small internal reversal to break symmetry
        if (idx + round_id) % 2 == 1:
          block = list(reversed(block))
        blocks.append(block)

      # interleave blocks; alternate order by round to increase mixing
      micro = []
      maxlen = max(len(b) for b in blocks) if blocks else 0
      order = list(range(len(blocks)))
      if round_id % 2 == 1:
        order.reverse()
      for i in range(maxlen):
        for j in order:
          blk = blocks[j]
          if i < len(blk):
            micro.append(blk[i])

      # fractional-span connectors (localized) to tie micro windows together
      micro_connectors = [
        (glo + int(round(0.08 * G)), glo + int(round(0.30 * G))),
        (glo + int(round(0.60 * G)), glo + int(round(0.92 * G))),
        (glo + int(round(0.26 * G)), glo + int(round(0.56 * G))),
        (glo + int(round(0.44 * G)), glo + int(round(0.78 * G))),
      ]
      if alt:
        # add one longer connector in alternate phase
        micro_connectors.append((glo + int(round(0.18 * G)), glo + int(round(0.84 * G))))

      for a, b in micro_connectors:
        if b > a:
          micro.append((a, b))

      # trim to budget
      if len(micro) > budget:
        micro = micro[:budget]
      return micro

    # Execute primary micro round (conservative budget)
    remaining = CAP - len(T)
    if remaining > 12:
      bud = max(8, remaining // 3)
      mic1 = micro_windows_round(T, round_id=0, budget=bud, alt=False)
      if mic1:
        take = min(len(mic1), CAP - len(T))
        T.extend(mic1[:take])

    # Execute alternate micro round (uses remaining room)
    remaining = CAP - len(T)
    if remaining > 8:
      mic2 = micro_windows_round(T, round_id=1, budget=min(128, remaining), alt=True)
      if mic2:
        take = min(len(mic2), CAP - len(T))
        T.extend(mic2[:take])

  # Final small set of long-range cross-scale connectors to force late-color creation
  if len(T) < CAP:
    lo = min(l for l, r in T)
    hi = max(r for l, r in T)
    span = max(1, hi - lo)
    final_cross = [
      (lo + int(0.02 * span), hi - int(0.02 * span)),
      (lo + int(0.12 * span), hi - int(0.12 * span)),
      (lo + int(0.28 * span), hi - int(0.28 * span)),
    ]
    for c in final_cross:
      if len(T) >= CAP:
        break
      if c[1] > c[0]:
        T.append(c)

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