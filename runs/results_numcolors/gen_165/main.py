# EVOLVE-BLOCK-START

def construct_intervals():
  """
  Deterministic round-robin parity spine with two thin micro-phases.
  Returns a list of (l, r) integer open intervals for FirstFit.
  """
  # Tight capacity to stay near the empirically best regime; keeps < 10000
  CAP = 9600

  # Fixed four-template bank (round-robin)
  TEMPLATES = [
    (2, 6, 10, 14),  # T1: classic KT
    (1, 5, 9, 13),   # T2: left-shifted
    (3, 7, 11, 15),  # T3: right-shifted
    (4, 8, 12, 16),  # T4: alternate wide
  ]

  # Seed: one unit interval; six rounds fit comfortably under CAP in this scheme.
  intervals = [(0, 1)]

  def spine_round(curr, starts, ridx, add_cross4=True, interleave_even=True):
    """
    Apply one parity-aware round using the given start template.
    connectors: classic + optional cross4.
    Even rounds: interleave blocks; odd rounds: sequential.
    """
    lo = min(l for l, _ in curr)
    hi = max(r for _, r in curr)
    delta = max(1, hi - lo)

    # Build four translated blocks
    blocks = []
    for s in starts:
      base = s * delta - lo
      blk = [(l + base, r + base) for (l, r) in curr]
      blocks.append(blk)

    # Interleave on even rounds; sequential on odd rounds
    out = []
    if interleave_even and (ridx % 2 == 0):
      m = max(len(b) for b in blocks)
      for i in range(m):
        for b in blocks:
          if i < len(b):
            out.append(b[i])
    else:
      for b in blocks:
        out.extend(b)

    # Classic four connectors
    s0, s1, s2, s3 = starts
    connectors = [
      ((s0 - 1) * delta, (s1 - 1) * delta),  # left cap
      ((s2 + 2) * delta, (s3 + 2) * delta),  # right cap
      ((s0 + 2) * delta, (s2 - 1) * delta),  # cross 1
      ((s1 + 2) * delta, (s3 - 1) * delta),  # cross 2
    ]
    # Add a single cross4 connector on odd rounds only (sparser -> safer for omega)
    if add_cross4 and (ridx % 2 == 1):
      connectors.append(((s0 + 4) * delta, (s3 + 4) * delta))

    out.extend(connectors)
    return out

  # Stage 1: six spine rounds, round-robin templates with parity behavior.
  for ridx in range(6):
    predicted = 4 * len(intervals) + 4
    if predicted > CAP:
      break
    starts = TEMPLATES[ridx % 4]
    intervals = spine_round(intervals, starts, ridx, add_cross4=True, interleave_even=True)

  # If already near the cap, return the backbone.
  if len(intervals) >= CAP - 32:
    return intervals

  # Stage 2: Two ultra-thin micro-phases at divisors 2 and 3.
  # Thin sampling + parity interleaving; minimal connectors only.
  def thin_seed(curr, k):
    n = len(curr)
    if n == 0 or k <= 0:
      return []
    step = max(1, n // k)
    return curr[::step][:k]

  def micro_phase(curr, rid, scale_div, max_budget):
    if not curr or max_budget <= 0:
      return []
    glo = min(l for l, _ in curr)
    ghi = max(r for _, r in curr)
    G = max(1, ghi - glo)
    delta = max(1, G // scale_div)

    # Choose a very thin seed to limit growth and omega impact
    per_block = max(12, min(48, max_budget // 16))
    U = thin_seed(curr, per_block)
    if not U:
      return []

    starts = TEMPLATES[0]  # lock to classic for micro to stabilize geometry
    ulo = min(l for l, _ in U)

    # Build four translated blocks with internal parity reversal
    blocks = []
    for s in starts:
      base = s * delta - ulo
      blk = [(l + base, r + base) for (l, r) in U]
      if ((s // 2) % 2) == (rid % 2):
        blk = list(reversed(blk))
      blocks.append(blk)

    out = []
    # Parity interleaving at micro-level
    order = blocks if (rid % 2 == 0) else list(reversed(blocks))
    m = max(len(b) for b in blocks)
    for i in range(m):
      for b in order:
        if i < len(b):
          out.append(b[i])

    # Minimal connectors (localized)
    s0, s1, s2, s3 = starts
    connectors = [
      ((s0 - 1) * delta + glo, (s1 - 1) * delta + glo),
      ((s2 + 2) * delta + glo, (s3 + 2) * delta + glo),
      ((s0 + 2) * delta + glo, (s2 - 1) * delta + glo),
    ]
    # Gate one cross4 in micro only on rid == 1 (odd parity), to stay sparse
    if rid % 2 == 1:
      connectors.append(((s0 + 4) * delta + glo, (s3 + 4) * delta + glo))
    out.extend(connectors)

    # Enforce budget
    if len(out) > max_budget:
      out = out[:max_budget]
    return out

  # Apply two micro-phases with tight budgets
  remaining = CAP - len(intervals)
  if remaining > 0:
    mp1 = micro_phase(intervals, rid=0, scale_div=2, max_budget=remaining // 2)
    if mp1:
      add = mp1[:CAP - len(intervals)]
      intervals.extend(add)

  remaining = CAP - len(intervals)
  if remaining > 0:
    mp2 = micro_phase(intervals, rid=1, scale_div=3, max_budget=remaining)
    if mp2:
      add = mp2[:CAP - len(intervals)]
      intervals.extend(add)

  return intervals

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()