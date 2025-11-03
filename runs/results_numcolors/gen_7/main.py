# EVOLVE-BLOCK-START

def construct_intervals(rounds=3):
  """
  Improved recursive wave construction inspired by arXiv:1506.00192.

  Key ideas:
  - Alternate clone order (reverse every other block) to disrupt FirstFit reuse.
  - Inject a few small "preload" blocker intervals at the start of each block
    so low colors become occupied early and clones are forced to seek higher colors.
  - Add short bridges between adjacent blocks to couple their color choices.
  - Use sparse caps and occasional long sentries to reserve colors across rounds.
  - Cap the number of rounds to avoid explosive growth.

  These arrival-order engineering techniques aim to increase FirstFit's color
  usage while keeping the offline optimum (largest clique) small by spacing blocks
  and limiting how many caps overlap a single point.
  """
  T = [(0.0, 1.0)]
  # safety cap on rounds to avoid creating enormous sequences
  rounds = max(1, min(rounds, 4))
  for round_idx in range(rounds):
    lo = min(l for l, r in T)
    hi = max(r for l, r in T)
    delta = hi - lo
    S = []
    # Widely spaced block starts; add more blocks in later rounds
    starts = [2, 6, 10, 14]
    if round_idx >= 1:
      starts += [18, 22, 26, 30]
    if round_idx >= 2:
      starts += [34, 38]
    # For each block: prepend a few short "preload" blockers, then append a cloned
    # copy of T (alternating order). Add a short bridge to the next block to couple colors.
    preload_blockers = 3
    for idx, start in enumerate(starts):
      for b in range(preload_blockers):
        off = 0.08 * b
        # place blocker inside the block interval [delta*start, delta*(start+1)]
        S.append((delta * (start + 0.08 + off), delta * (start + 0.08 + off + 0.18)))
      base = T if (idx % 2 == 0) else list(reversed(T))
      for l, r in base:
        S.append((delta * start + l - lo, delta * start + r - lo))
      # short bridge to next block (only overlaps two neighbouring blocks)
      if idx + 1 < len(starts):
        nxt = starts[idx + 1]
        S.append((delta * (start + 0.6), delta * (nxt + 0.35)))
    # local gadget: keeps omega small but increases FF pressure
    S += [
      (delta * 1, delta * 5),
      (delta * 12, delta * 16),
      (delta * 4, delta * 9),
      (delta * 8, delta * 13),
    ]
    # sparse caps (moderate spans), spaced to avoid creating large cliques
    max_start = max(starts)
    for j in range(5, int(max_start), 6):
      S.append((delta * (j - 1.5), delta * (j + 1.5)))
    # occasional long sentries to occupy low colors early (only in first round)
    if round_idx == 0:
      S += [(delta * 0.5, delta * 1.5), (delta * 1.2, delta * 2.2)]
    T = S
  return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()