# EVOLVE-BLOCK-START

def construct_intervals(rounds=3, seed_lo=0.0, seed_scale=1.0):
  """
  Deterministic, adaptive multi-round wave construction.
  - rounds: number of recursive expansion rounds
  - seed_lo / seed_scale: deterministic knobs for initial placement and scaling
  Returns:
    intervals: list of tuples (l, r) representing open intervals
  """
  # Start with a compact core
  T = [(seed_lo, seed_lo + 1.0 * seed_scale)]
  for round_idx in range(rounds):
    # Furthest-left/right boundaries of current wave
    lo = min(l for l, r in T)
    hi = max(r for l, r in T)
    delta = hi - lo if hi > lo else 1.0  # avoid degenerate

    S = []
    # Adaptive starts: denser in early rounds, controlled growth later
    if round_idx == 0:
      starts = [2, 4, 6, 8, 10, 12, 14, 16]
    else:
      # add a few extra starts but cap total to prevent exponential blowup
      starts = [2, 4, 6, 8, 10, 12, 14, 16, 18, 22, 26]

    # Create translated copies of current T
    for i, start in enumerate(starts):
      # Alternate inner order to reduce color re-use
      base = T if (i % 2 == 0) else list(reversed(T))
      S += [(delta * start + l - lo, delta * start + r - lo) for l, r in base]

      # Optional light bridging to couple colors across adjacent blocks
      if i + 1 < len(starts):
        nxt = starts[i + 1]
        # bridge spans inside current to inside next block
        S.append((delta * (start + 0.5), delta * (nxt + 0.5)))

    # Local four-interval gadget to keep omega small but push FF usage
    S += [
      (delta * 1.0, delta * 5.0),
      (delta * 12.0, delta * 16.0),
      (delta * 4.0, delta * 9.0),
      (delta * 8.0, delta * 13.0)
    ]

    # Tiny caps to further push color usage without large omega
    caps = []
    cap1 = max(1.0, delta * 0.4)
    cap2 = delta * 1.8
    # place a few small caps at spaced positions
    for j in range(3):
      caps.append((cap1 * (j + 1), cap1 * (j + 1) + cap2 * 0.2))
    S += caps

    T = S

  return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()