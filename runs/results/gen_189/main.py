# EVOLVE-BLOCK-START

def construct_intervals(iterations=4, extra_first=True):
  """
  Construct a sequence of intervals for the FirstFit adversary.

  Improvements over the canonical variant:
    - cycle among four offset patterns (A,B,C,D) instead of only two,
      creating more varied tilings across levels;
    - apply a tiny deterministic perturbation to each copy's placement
      (depends on level and start) to break perfect symmetry and induce
      additional cross-level overlaps;
    - on alternating levels insert the four long blockers interleaved
      among the copies (and reverse the copy order on those levels)
      so that blockers appear earlier relative to some copies and
      pollute small colors more effectively;
    - deduplicate identical intervals (preserve arrival order) to keep
      the instance compact;
    - finally normalize endpoints to a compact even-integer grid.

  These are conservative, deterministic changes intended to increase
  the number of colors FirstFit uses while keeping omega controlled.
  """

  # four cyclic offset patterns (A,B,C,D)
  offsets_patterns = [
    (2, 6, 10, 14),   # A: canonical
    (1, 5, 9, 13),    # B: shifted
    (3, 7, 11, 15),   # C: shifted-right
    (0, 4, 8, 12)     # D: leftmost
  ]

  # canonical blockers (kept same geometric multipliers)
  blockers_template = [(1, 5), (12, 16), (4, 9), (8, 13)]

  # start with a tiny seed to keep clique small
  T = [(0.0, 1.0)]

  for lvl in range(iterations):
    lo = min(l for l, r in T)
    hi = max(r for l, r in T)
    delta = hi - lo
    S = []

    # choose base offsets from the 4-pattern cycle
    base_offsets = offsets_patterns[lvl % len(offsets_patterns)]
    if extra_first and lvl == 0:
      offs = tuple(list(base_offsets) + [18])  # small extra coupling on first level
    else:
      offs = base_offsets

    # decide whether to interleave blockers on this level (alternate levels)
    interleave = (lvl % 2 == 1)
    # when interleaving, also reverse copy order to vary arrival
    reverse_order = interleave

    ordered_offs = list(offs)
    if reverse_order:
      ordered_offs = list(reversed(ordered_offs))

    # small deterministic perturbation generator (range ≈ [-0.12, +0.12])
    def perturb_for(start):
      # integerize start to keep perturb deterministic and small
      key = int(start)
      v = ((lvl * 3 + key) % 5) - 2  # in {-2,-1,0,1,2}
      return v * 0.06

    # Build S by adding copies per start and optionally interleaving a blocker
    for idx, start in enumerate(ordered_offs):
      p = perturb_for(start)
      offset = delta * (start + p) - lo
      # append translated copy of entire previous level T for this start
      for (l, r) in T:
        S.append((offset + l, offset + r))
      # interleave the corresponding blocker (if any) to affect arrival order
      if interleave and idx < len(blockers_template):
        a, b = blockers_template[idx]
        S.append((delta * a, delta * b))

    # if not interleaving, append blockers after all copies (canonical order)
    if not interleave:
      for (a, b) in blockers_template:
        S.append((delta * a, delta * b))

    T = S

  # deduplicate identical intervals while preserving arrival order (keeps instance compact)
  seen = set()
  T_dedup = []
  for (l, r) in T:
    key = (float(l), float(r))
    if key not in seen:
      seen.add(key)
      T_dedup.append((l, r))
  T = T_dedup

  # Normalize endpoints to a compact integer grid while preserving order.
  # Map each unique endpoint to an increasing even integer (to ensure
  # positive lengths and avoid degeneracy).
  endpoints = sorted(set([x for seg in T for x in seg]))
  coord = {}
  cur = 0
  for e in endpoints:
    coord[e] = cur
    cur += 2
  normalized = [(coord[l], coord[r]) for (l, r) in T]

  return normalized

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()