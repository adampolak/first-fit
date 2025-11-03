# EVOLVE-BLOCK-START

def construct_intervals(iterations: int = 4):
  """
  Adversarial interval sequence construction:
    - Build a 4-copy-style recursive backbone with cycling offsets and tiny deterministic perturbations.
    - Normalize the backbone and then deterministically append "waves" that increase FirstFit
      without increasing the clique number (OPT).
    - Aggressively prune intervals that do not affect FirstFit or OPT.
    - Return endpoints remapped to an increasing even integer grid.

  Arguments:
    iterations: recursion depth hint (default 4)
  Returns:
    list of (l, r) tuples (open intervals), integer endpoints
  """

  # ---------- Utilities ----------
  from math import floor

  # small deterministic "pseudo-random" fraction generator (no randomness)
  def det_frac(a, b):
    # returns fraction in (-0.25, 0.25) deterministically
    return (((a * 1315423911) ^ (b * 2654435761)) % 1000 - 500) / 2000.0 * 0.25

  # Convert list of intervals into a hashable key
  def _keyify(intervals):
    return tuple((int(l), int(r)) if (isinstance(l, (int,)) and isinstance(r, (int,))) else (float(l), float(r)) for (l, r) in intervals)

  # Memoized clique and FirstFit
  _clique_cache = {}
  _ff_cache = {}

  def clique_number(intervals):
    """Sweep algorithm for open intervals. Right endpoints processed before left at ties."""
    key = _keyify(intervals)
    if key in _clique_cache:
      return _clique_cache[key]
    events = []
    for (l, r) in intervals:
      if l < r:
        events.append((l, +1))
        events.append((r, -1))
    events.sort(key=lambda e: (e[0], 0 if e[1] == -1 else 1))
    cur = best = 0
    for _, t in events:
      cur += t
      if cur > best:
        best = cur
    _clique_cache[key] = best
    return best

  def firstfit_colors(intervals):
    """Simulate FirstFit on arrival order using per-color last end tracking."""
    key = _keyify(intervals)
    if key in _ff_cache:
      return _ff_cache[key]
    last_end = []  # right endpoint of last interval in each color
    for (l, r) in intervals:
      placed = False
      for i, le in enumerate(last_end):
        # open intervals -> no overlap if l >= le
        if l >= le:
          last_end[i] = r
          placed = True
          break
      if not placed:
        last_end.append(r)
    res = len(last_end)
    _ff_cache[key] = res
    return res

  def normalize_grid(intervals):
    """Map all unique endpoints to compact even integers (0,2,4,...)."""
    if not intervals:
      return []
    pts = sorted(set([x for seg in intervals for x in seg]))
    mapping = {v: 2 * i for i, v in enumerate(pts)}
    return [(mapping[l], mapping[r]) for (l, r) in intervals]

  # ---------- Core backbone builder ----------
  def build_backbone(depth, cycle_patterns, translation_mode='left'):
    """
    Build a recursive backbone using a cycle of offset patterns.
    - cycle_patterns: list/tuple of offset-tuples to cycle among levels
    - translation_mode: 'left' or 'center' anchoring (alternates can be simulated by caller)
    """
    T = [(0.0, 1.0)]
    for lvl in range(depth):
      lo = min(l for l, r in T)
      hi = max(r for l, r in T)
      delta = hi - lo
      center = (lo + hi) / 2.0

      pattern = cycle_patterns[lvl % len(cycle_patterns)]
      # Deterministic tiny perturbation per (lvl, copy_index)
      copies = []
      for ci, base in enumerate(pattern):
        # small deterministic fractional shift to break symmetry
        pf = det_frac(lvl + 1, ci + 1) * 0.15  # up to ~0.0375 * base scale
        start = base + pf
        if translation_mode == 'left':
          off = delta * start - lo
        else:
          off = delta * start - center
        # Alternate reverse ordering of the inner sequence to vary arrival order
        seq = list(T) if ((lvl + ci) % 2 == 0) else list(reversed(T))
        for (l, r) in seq:
          copies.append((l + off, r + off))

      # Four canonical blockers (Figure-4 style), anchored left relative to current window
      blockers = [
        (delta * 1.0,  delta * 5.0),
        (delta * 12.0, delta * 16.0),
        (delta * 4.0,  delta * 9.0),
        (delta * 8.0,  delta * 13.0),
      ]
      # Alternate schedule deterministically: after vs split on even/odd lvl
      if lvl % 2 == 0:
        S = copies + blockers
      else:
        h = max(1, len(pattern) // 2)
        # split roughly: first h copies, blockers, remaining copies
        first_block = copies[:h * len(T)]
        second_block = copies[h * len(T):]
        S = first_block + blockers + second_block
      # Keep floats here; normalized later.
      T = S
    return T

  # ---------- Wave injection ----------
  def inject_waves(base_intervals, max_waves=200, try_centers=80):
    """
    Deterministically attempt to append waves. Each candidate wave is accepted only if:
      - it does not increase the clique number (OPT stays the same)
      - it strictly increases the FirstFit color count
    Waves are appended to the end (arrive after backbone).
    """
    # Work on a mutable copy (intervals are integerized to keep search stable)
    current = list(base_intervals)
    om_base = clique_number(current)
    cols = firstfit_colors(current)
    if om_base <= 0:
      return current, cols, om_base

    lo = min(l for l, r in current)
    hi = max(r for l, r in current)
    rng = hi - lo
    if rng <= 0:
      return current, cols, om_base

    # Choose candidate centers at deterministic spacing
    spacing = max(2, (rng // max(4, try_centers)))
    centers = []
    start = lo + spacing // 2
    pos = start
    while pos < hi - 1 and len(centers) < try_centers:
      centers.append(pos)
      pos += spacing

    # Wave lengths: odd small integers up to a fraction of the full span.
    max_len = max(1, rng // 6)
    wlens = [w for w in range(1, max(2, max_len + 1), 2)]
    # Try longer odd lengths last (so we prefer short waves first)
    wlens = sorted(wlens)

    inserted = 0
    # Greedy deterministic scan
    for c in centers:
      if inserted >= max_waves:
        break
      accepted = False
      for w in wlens:
        lw = c - w // 2
        rw = lw + w
        if lw >= rw:
          continue
        # Small safety: keep wave inside bounding box with 1 unit padding
        if lw < lo - 1 or rw > hi + 1:
          continue
        cand = current + [(lw, rw)]
        # Evaluate candidate
        new_om = clique_number(cand)
        if new_om > om_base:
          # would increase OPT; forbidden
          continue
        new_cols = firstfit_colors(cand)
        if new_cols > cols:
          # Accept wave: it forces a new color without increasing OPT
          current = cand
          cols = new_cols
          inserted += 1
          accepted = True
          # small optimization: update caches may be beneficial; they will be updated by functions
          break
      # proceed to next center
      # (We don't accept waves that don't immediately increase FF; that keeps control of OPT)
    return current, cols, om_base

  # ---------- Deterministic pruning ----------
  def prune_preserving_metrics(intervals):
    """
    Repeatedly remove single intervals (left-to-right) whenever both
    clique number and FirstFit color count remain unchanged.
    Deterministic greedy loop until no change.
    """
    T = list(intervals)
    if not T:
      return T
    om = clique_number(T)
    cols = firstfit_colors(T)
    changed = True
    # To keep runtime bounded, limit total removals to a fraction of size
    max_removals = max(1, len(T) // 2)
    removals = 0
    while changed and removals < max_removals:
      changed = False
      # iterate left-to-right to prefer removing later helpers
      i = 0
      while i < len(T):
        cand = T[:i] + T[i+1:]
        if not cand:
          i += 1
          continue
        if clique_number(cand) == om and firstfit_colors(cand) == cols:
          T = cand
          changed = True
          removals += 1
          # do not increment i: we've shifted list left
        else:
          i += 1
        if removals >= max_removals:
          break
    return T

  # ---------- Candidate exploration ----------
  # offsets patterns to cycle among levels (backbone-cycling)
  P_A = (2.0, 6.0, 10.0, 14.0)
  P_B = (1.0, 5.0, 9.0, 13.0)
  P_C = (3.0, 7.0, 11.0, 15.0)
  P_D = (0.0, 4.0, 8.0, 12.0)
  cycle_options = [
    (P_A, P_B, P_C, P_D),
    (P_A, P_B),
    (P_A, P_C),
  ]
  translation_modes = ['left', 'center']
  depth_candidates = [max(2, iterations - 1), iterations, min(iterations + 1, 5)]
  depth_candidates = sorted(set([d for d in depth_candidates if 2 <= d <= 5]))

  best = None  # (ratio, -n, cols, om, intervals)

  # Limit search size modestly to keep runtime reasonable
  for cycle in cycle_options:
    for trans in translation_modes:
      for depth in depth_candidates:
        # Build backbone
        core = build_backbone(depth, cycle, translation_mode=trans)
        # Normalize the backbone before wave injection for stability
        core_norm = normalize_grid(core)
        # Evaluate baseline
        om_core = clique_number(core_norm)
        if om_core <= 0:
          continue
        cols_core = firstfit_colors(core_norm)
        # Inject waves (deterministically)
        candidate_with_waves, cols_after, om_after = inject_waves(core_norm, max_waves=250, try_centers=96)
        # If injection did not change OM (it should not), compute metrics
        final_om = om_core  # we kept waves only when OM unchanged
        final_cols = cols_after
        final_intervals = list(candidate_with_waves)
        # Aggressive deterministic pruning
        pruned = prune_preserving_metrics(final_intervals)
        pruned = normalize_grid(pruned)  # compact final representation
        final_om = clique_number(pruned)
        final_cols = firstfit_colors(pruned)
        n_final = len(pruned)
        if final_om <= 0:
          continue
        ratio = final_cols / final_om
        cand = (ratio, -n_final, final_cols, final_om, pruned)
        if best is None:
          best = cand
        else:
          # prefer larger ratio; tie-break by fewer intervals, then more colors
          if cand[0] > best[0] + 1e-12 or (abs(cand[0] - best[0]) <= 1e-12 and (cand[1] > best[1] or (cand[1] == best[1] and cand[2] > best[2]))):
            best = cand

  # Fallback: canonical Figure-4 4-iteration baseline (safe)
  if best is None:
    T = [(0.0, 1.0)]
    for _ in range(max(3, iterations)):
      lo = min(l for l, r in T)
      hi = max(r for l, r in T)
      delta = hi - lo
      S = []
      for start in (2.0, 6.0, 10.0, 14.0):
        off = delta * start - lo
        for (l, r) in T:
          S.append((l + off, r + off))
      S += [
        (delta * 1.0,  delta * 5.0),
        (delta * 12.0, delta * 16.0),
        (delta * 4.0,  delta * 9.0),
        (delta * 8.0,  delta * 13.0),
      ]
      T = S
    return normalize_grid(T)

  # Return best found (already normalized)
  return best[4]

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()