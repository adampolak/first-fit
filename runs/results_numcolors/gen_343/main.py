# EVOLVE-BLOCK-START

def construct_intervals(enable_alt_microphase=True):
  """
  Greedy, simulator-driven interval construction designed to maximize
  FirstFit colors while keeping the offline clique (omega) <= 10.

  Args:
    enable_alt_microphase (bool): preserved parameter for compatibility.

  Returns:
    List of integer open intervals [(l,r), ...] in the online arrival order.
  """

  # Hard capacity guard
  CAP = 9800

  # Hard offline optimum bound we enforce
  OMEGA_LIMIT = 10

  # Deterministic "random" helper using a simple LCG so results are reproducible
  class DeterministicRNG:
    def __init__(self, seed=1234567):
      self.x = seed
    def rand(self):
      # 32-bit LCG
      self.x = (1103515245 * self.x + 12345) & 0x7fffffff
      return self.x
    def randint(self, a, b):
      return a + (self.rand() % (b - a + 1))

  rng = DeterministicRNG(seed=314159265)

  # A compact integer grid representing the real line positions we'll use.
  # We will keep endpoints integer, and test coverage only at integer sample points.
  # We make the span moderate so we can place many intervals without huge coordinates.
  GRID_SPAN = 20000  # total coordinate span allowed (will shift/scatter inside it)
  # We'll partition the grid into "columns" and "mixers" to get structure.
  # Keep a global offset to avoid negative endpoints.
  OFFSET = 1000

  # Internal representation: list of (l, r) integer tuples
  intervals = []

  # Fast coverage structure for integer points using dictionary: count per integer point
  coverage = {}

  def inc_coverage_for_interval(iv):
    l, r = iv
    # we treat open intervals (l, r) but for integer grid we increment for p in [l, r-1]
    for p in range(l, r):
      coverage[p] = coverage.get(p, 0) + 1

  def dec_coverage_for_interval(iv):
    l, r = iv
    for p in range(l, r):
      c = coverage.get(p, 0)
      if c <= 1:
        # remove key if becomes 0
        if p in coverage:
          del coverage[p]
      else:
        coverage[p] = c - 1

  # Simulate FirstFit coloring for a sequence of intervals (online).
  # Returns list colors[0..n-1] where colors[i] is color assigned (1-based).
  def simulate_firstfit(seq):
    # For speed, we construct an interval list and perform naive online FF.
    # We'll maintain for each active color the rightmost endpoint (so we know if color available).
    # But intervals may overlap non-nestedly; simplest safe implementation:
    colors = []
    # maintain current intervals' endpoints by color to test conflicts
    color_intervals = {}  # color -> list of active intervals (l,r) that may conflict with future intervals
    # For correctness we simply check overlaps with already assigned intervals of each color
    for (l, r) in seq:
      c = 1
      while True:
        # check whether color c is usable (i.e., it has no assigned interval that intersects (l,r))
        conflict = False
        if c in color_intervals:
          for (al, ar) in color_intervals[c]:
            # intervals open: overlap iff not (ar <= l or r <= al)
            if not (ar <= l or r <= al):
              conflict = True
              break
        if not conflict:
          # assign c
          colors.append(c)
          color_intervals.setdefault(c, []).append((l, r))
          break
        c += 1
    return colors

  # Helper: max coverage count across integer sample points (current omega observed)
  def current_omega():
    if not coverage:
      return 0
    return max(coverage.values())

  # Helper: feasibility check: would adding iv violate OMEGA_LIMIT at any integer point?
  def feasible_add(iv):
    l, r = iv
    for p in range(l, r):
      if coverage.get(p, 0) + 1 > OMEGA_LIMIT:
        return False
    return True

  # Helper: normalize interval endpoints (ensure r > l and inside allowed GRID_SPAN)
  def make_interval(l, r):
    if r <= l:
      r = l + 1
    # clamp into span [OFFSET, OFFSET+GRID_SPAN)
    L = max(OFFSET, min(OFFSET + GRID_SPAN - 1, l))
    R = max(L + 1, min(OFFSET + GRID_SPAN, r))
    return (L, R)

  # Seed construction: build a small scaffold of "columns" (stacked intervals) to set omega baseline.
  # We place C columns, each with H stacked intervals (so clique at column points = H).
  H = 8  # base stack height (leave a little headroom up to 10 for micro phases)
  C = 40  # number of columns to create initially; tradeoff space vs capacity
  SPACING = 40  # spacing between columns
  col_positions = []
  base_x = OFFSET + 50
  for j in range(C):
    pos = base_x + j * SPACING
    col_positions.append(pos)

  # Create stacks: in each column, place H overlapping intervals that all include [pos, pos+1)
  for pos in col_positions:
    # stack intervals with slightly varying endpoints to ensure they pairwise overlap at pos
    for level in range(H):
      l = pos - (2 + level % 2)  # slight jitter by level
      r = pos + 5 + (level % 3)
      iv = make_interval(l, r)
      if feasible_add(iv):
        intervals.append(iv)
        inc_coverage_for_interval(iv)
      else:
        # Shouldn't typically happen at seed time given small counts
        pass
      if len(intervals) >= CAP:
        # enforce cap early
        return intervals

  # At this point offline omega is at most H (≈8). We'll try to force FirstFit to use many colors.
  # We'll maintain representatives: for each color value currently present in FF coloring of intervals,
  # we will keep at least one interval (a representative) located in a distinct column, so later
  # candidate intervals can touch all representatives without creating large cliques.
  # We'll keep representatives updated after each insertion by simulating FF on current sequence.
  def build_representative_map(seq, colors):
    # Map color -> index of representative interval in seq
    rep = {}
    for idx, c in enumerate(colors):
      if c not in rep:
        rep[c] = idx
    return rep

  # compute current FirstFit coloring
  ff_colors = simulate_firstfit(intervals)
  rep_map = build_representative_map(intervals, ff_colors)
  current_ff_max = max(ff_colors) if ff_colors else 0

  # Greedy simulator-driven phase:
  # We repeatedly attempt to find small candidate intervals that, when added, increase the FF max color.
  # Candidate generation is deterministic, explores many short intervals centered at column positions and
  # at inter-column midpoints. We limit search breadth per iteration to keep runtime bounded.
  ITER_LIMIT = CAP - len(intervals)
  # Parameters controlling search efforts
  CANDIDATES_PER_ITER = 350  # number of candidates to evaluate per greedy iteration
  MAX_TRIES_NO_IMPROVE = 220  # if we cannot find improvement for this many iterations, we fall back

  no_improve_count = 0
  iter_count = 0
  while len(intervals) < CAP and no_improve_count < MAX_TRIES_NO_IMPROVE and iter_count < ITER_LIMIT:
    iter_count += 1
    improved = False
    best_candidate = None
    best_new_max = current_ff_max
    # Deterministic sampling pattern: sample centers among column positions and inter-column midpoints,
    # plus a few random-looking offsets from rng to diversify.
    candidates = []
    # sample near columns
    for j in range(min(C, max(1, CANDIDATES_PER_ITER // 5))):
      idx = (iter_count + j) % len(col_positions)
      cen = col_positions[idx]
      # produce a few short lengths
      for length in (1, 2, 3, 4):
        l = cen - length
        r = cen + 1
        candidates.append(make_interval(l, r))
    # sample between columns
    for j in range(min(C - 1, max(1, CANDIDATES_PER_ITER // 5))):
      idx = (iter_count + 2 * j) % (len(col_positions) - 1)
      cen = (col_positions[idx] + col_positions[idx + 1]) // 2
      for length in (1, 2, 3):
        l = cen - 1
        r = cen + length
        candidates.append(make_interval(l, r))
    # sample some increasingly wide intervals that bridge many columns (but short enough to keep clique)
    for wid in (6, 10, 14, 18):
      # pick deterministic start shifting
      start = OFFSET + (iter_count * 17) % (GRID_SPAN - wid - 10)
      candidates.append(make_interval(start, start + wid))
    # deterministic pseudo-random candidates (small number)
    for t in range(max(1, CANDIDATES_PER_ITER - len(candidates))):
      a = rng.randint(OFFSET + 2, OFFSET + GRID_SPAN - 10)
      b = a + rng.randint(1, 6)
      candidates.append(make_interval(a, b))

    # Evaluate candidates deterministically; keep the one that maximizes FF max color upon insertion
    # but respect OMEGA_LIMIT.
    base_seq = intervals  # current sequence
    base_colors = ff_colors
    # We will early prune: if candidate overlaps any integer point exceeding OMEGA_LIMIT, skip.
    for iv in candidates:
      if not feasible_add(iv):
        continue
      # Quick simulate: append iv to sequence and run FF simulation by incremental method:
      # We simulate full FirstFit on base_seq + [iv].
      # For speed we call simulate_firstfit on a short copy.
      test_seq = base_seq + [iv]
      test_colors = simulate_firstfit(test_seq)
      test_max = max(test_colors) if test_colors else 0
      if test_max > best_new_max:
        best_new_max = test_max
        best_candidate = iv
        # We can break early if we make substantial progress (greedy)
        # But keep scanning some more to find possibly better immediate jump.
    if best_candidate is not None:
      # Accept best candidate
      intervals.append(best_candidate)
      inc_coverage_for_interval(best_candidate)
      ff_colors = simulate_firstfit(intervals)
      rep_map = build_representative_map(intervals, ff_colors)
      current_ff_max = max(ff_colors) if ff_colors else 0
      improved = True
      no_improve_count = 0
    else:
      no_improve_count += 1
      # If no candidate improved FF, we add a controlled long connector to mix colors deterministically.
      # The connector spans a moderate span capturing many representatives but arranged so its insertion
      # does not break OMEGA_LIMIT. We build connectors between selected columns spaced apart.
      # Connector length grows slowly with iter_count to diversify interactions.
      # Choose leftmost column index depending on iter_count
      left_idx = (iter_count * 3) % max(1, len(col_positions) - 6)
      right_idx = min(len(col_positions) - 1, left_idx + 6 + (iter_count % 4))
      l = col_positions[left_idx] - 3
      r = col_positions[right_idx] + 6 + (iter_count % 5)
      conn = make_interval(l, r)
      if feasible_add(conn):
        intervals.append(conn)
        inc_coverage_for_interval(conn)
        ff_colors = simulate_firstfit(intervals)
        rep_map = build_representative_map(intervals, ff_colors)
        current_ff_max = max(ff_colors) if ff_colors else 0
        improved = True
        no_improve_count = 0
      else:
        # as last resort, add a very short pin in the least covered area to expand search base
        # find a point p with smallest coverage
        # deterministic scan across sample points to find p with minimal coverage
        min_cov = OMEGA_LIMIT
        min_p = OFFSET + 2
        # sample across a fraction of grid for speed
        stride = max(1, GRID_SPAN // 400)
        for p in range(OFFSET + 1, OFFSET + GRID_SPAN - 1, stride):
          c = coverage.get(p, 0)
          if c < min_cov:
            min_cov = c
            min_p = p
            if min_cov == 0:
              break
        pin = make_interval(min_p, min_p + 1)
        if feasible_add(pin):
          intervals.append(pin)
          inc_coverage_for_interval(pin)
          ff_colors = simulate_firstfit(intervals)
          rep_map = build_representative_map(intervals, ff_colors)
          current_ff_max = max(ff_colors) if ff_colors else 0
          improved = True
          no_improve_count = 0
        else:
          # cannot add anything reasonable: break out
          break

    # continue loop while we still have capacity
    if not improved:
      # safety guard: avoid infinite loop
      break

  # After greedy phase, if still room, add several deterministic "mixers" and final pins to increase pressure further.
  while len(intervals) < CAP:
    # Add deterministic connectors spanning different windows but respecting omega
    base_lo = OFFSET + 20
    shift = (len(intervals) * 37) % (GRID_SPAN // 3)
    L = base_lo + shift
    R = L + 50 + (len(intervals) % 37)
    iv = make_interval(L, R)
    if feasible_add(iv):
      intervals.append(iv)
      inc_coverage_for_interval(iv)
    else:
      # add a tiny pin somewhere low-covered
      # deterministic scan for a low covered point in a small sample region
      got = False
      for p in range(OFFSET + 1, OFFSET + GRID_SPAN - 1, max(1, GRID_SPAN // 600)):
        if coverage.get(p, 0) < OMEGA_LIMIT:
          pin = make_interval(p, p + 1)
          if feasible_add(pin):
            intervals.append(pin)
            inc_coverage_for_interval(pin)
            got = True
            break
      if not got:
        break

  # Final normalization/trimming to CAP and ensure integer endpoints and positive lengths
  if len(intervals) > CAP:
    intervals = intervals[:CAP]

  normed = []
  for (l, r) in intervals:
    li = int(round(l))
    ri = int(round(r))
    if ri <= li:
      ri = li + 1
    normed.append((li, ri))

  return normed

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()