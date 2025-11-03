# EVOLVE-BLOCK-START

import math
from bisect import bisect_left, bisect_right

def _sweep_omega(intervals):
  """Compute maximum number of intervals covering any point (omega), using sweep-line."""
  if not intervals:
    return 0
  events = []
  for l, r in intervals:
    # open intervals (l,r): treat start as +1 at l, end as -1 at r
    events.append((l, 1))
    events.append((r, -1))
  events.sort()
  cur = 0
  best = 0
  for _, delta in events:
    cur += delta
    if cur > best:
      best = cur
  return best

def _firstfit_colors(intervals):
  """
  Simulate FirstFit on given sequence of intervals.
  Returns list of colors (1-based) corresponding to intervals.
  """
  colors = []
  # maintain for each color, a list of active intervals' indices or endpoints
  # naive approach: for each new interval check overlaps with all previous intervals.
  for i, (l, r) in enumerate(intervals):
    used = set()
    for j, (lj, rj) in enumerate(intervals[:i]):
      if lj < r and l < rj:
        used.add(colors[j])
    # find smallest positive integer not in used
    c = 1
    while c in used:
      c += 1
    colors.append(c)
  return colors

def _assigned_color_if_appended(intervals, colors, candidate):
  """
  Given current intervals and their assigned FirstFit colors, compute
  the color FirstFit would assign to candidate if appended next.
  This avoids recomputing all colors.
  """
  l, r = candidate
  used = set()
  for (li, ri), ci in zip(intervals, colors):
    if li < r and l < ri:
      used.add(ci)
  c = 1
  while c in used:
    c += 1
  return c

def _omega_if_appended(intervals, candidate):
  """Compute omega after appending candidate by a sweep merging two events arrays."""
  # faster to do a small sweep merging existing events with candidate events
  events = []
  for l, r in intervals:
    events.append((l, 1))
    events.append((r, -1))
  events.append((candidate[0], 1))
  events.append((candidate[1], -1))
  events.sort()
  cur = 0
  best = 0
  for _, delta in events:
    cur += delta
    if cur > best:
      best = cur
  return best

def _interesting_endpoints(intervals):
  """Return a small sorted list of distinct endpoints (deterministic sampling)."""
  if not intervals:
    return [0, 1]
  pts = sorted(set([x for l, r in intervals for x in (l, r)]))
  # sample endpoints: min, max, mid, quartiles (bounded to small list)
  n = len(pts)
  chosen = []
  idxs = [0, n-1, n//2, n//4, (3*n)//4]
  for idx in idxs:
    if 0 <= idx < n:
      chosen.append(pts[idx])
  # also include a handful of near-tail points for span candidates
  chosen.extend(pts[-6:])   # last few
  chosen.extend(pts[:6])    # first few
  # unique and sorted
  chosen = sorted(set(chosen))
  return chosen

def construct_intervals(seed_count=1, phase2_iters=2, enable_parabolic=False):
  """
  Adaptive, simulated-firstfit greedy construction.

  Inputs:
    - seed_count (int): ignored for complexity but preserved signature (keeps deterministic behavior).
    - phase2_iters (int): hint for number of "aggressive" greedy rounds (kept limited).
    - enable_parabolic (bool): kept for compatibility (unused in this method).

  Output:
    - list of (l, r) tuples (open intervals), in the order to be presented to FirstFit.
  """
  CAP = 9800
  MAX_OMEGA = 10  # keep offline optimum <= 10

  intervals = []

  # deterministic small seed to give structure
  # create a small symmetric anchor spine (deterministic)
  anchors = [(0, 1), (4, 5), (8, 9), (12, 13)]
  for a in anchors:
    if len(intervals) < CAP:
      intervals.append(a)

  # compute current FF colors
  colors = _firstfit_colors(intervals)

  # We will run a bounded number of greedy iterations to fill up to CAP.
  # Each iteration proposes a small deterministic candidate pool and picks the best
  # candidate (highest assigned color while preserving omega <= MAX_OMEGA).
  # The candidate pool is constructed from existing endpoints and a palette of lengths.
  LENS = [1, 2, 3, 5, 8, 13, 21, 34]  # Fibonacci-like sizes to mix granularities
  # Deterministic 'seed' parameter to make proposals repeatable
  iter_limit = min(CAP - len(intervals), 8200)  # keep an upper bound on iterations
  # phase2_iters increases the aggressiveness (more long spans at start)
  extra_long_factor = 1 + max(0, min(4, int(phase2_iters)))

  # Precompute existing events for faster omega checks occasionally
  # We'll recompute omega per candidate but we maintain a reasonable candidate pool size.
  for step in range(iter_limit):
    room = CAP - len(intervals)
    if room <= 0:
      break

    # current colors
    colors = _firstfit_colors(intervals)

    # gather interesting endpoints deterministically
    pts = _interesting_endpoints(intervals)

    candidates = []

    # 1) Short candidates anchored near interesting endpoints
    for p in pts:
      for L in LENS[:4]:
        l = p - L // 2
        r = l + L
        if r <= l:
          r = l + 1
        candidates.append((l, r))

    # 2) Mid candidates: span from one endpoint to another nearby one
    m = len(pts)
    for i in range(0, m, max(1, m // 8)):
      for j in range(i+1, min(m, i+6)):
        l = pts[i] - 1
        r = pts[j] + 1
        if r > l:
          candidates.append((l, r))

    # 3) Long candidates: occasional long spans, scaled by extra_long_factor
    span_lo = min(pts)
    span_hi = max(pts)
    span = max(1, span_hi - span_lo)
    for L in (int(span * 0.08), int(span * 0.15), int(span * 0.25)):
      if L <= 0:
        continue
      # derive a few anchored long placements
      for shift_frac in (0.02, 0.10, 0.28):
        l = span_lo + int(round(shift_frac * span)) - (L // 2)
        r = l + max(1, L) * extra_long_factor
        candidates.append((l, r))

    # 4) "Comb" candidates: span a few spaced anchors without being extremely long
    # use last few endpoints to craft varied spans
    tail = pts[-8:] if len(pts) >= 8 else pts
    for i in range(len(tail)):
      for j in range(i+1, min(len(tail), i+5)):
        l = tail[i] - 1
        r = tail[j] + 1
        candidates.append((l, r))

    # make candidate pool deterministic and unique
    seen = set()
    uniq_cands = []
    for c in candidates:
      # normalize to integers
      l, r = int(c[0]), int(c[1])
      if r <= l:
        r = l + 1
      key = (l, r)
      if key in seen:
        continue
      seen.add(key)
      uniq_cands.append(key)
      if len(uniq_cands) >= 120:  # cap candidate pool size for performance
        break

    # Evaluate candidates: compute assigned FF color and new omega
    best_score = (-1, 0, None)  # (assigned_color, -omega, interval)
    # To limit overhead, only fully evaluate limited candidates (top by heuristic length)
    # Heuristic: evaluate longer candidates first since they tend to hit many colors
    uniq_cands.sort(key=lambda iv: (iv[1] - iv[0]), reverse=True)

    # We'll evaluate up to N_EVAL candidates this round
    N_EVAL = 60
    for cand in uniq_cands[:N_EVAL]:
      assigned = _assigned_color_if_appended(intervals, colors, cand)
      new_omega = _omega_if_appended(intervals, cand)
      if new_omega > MAX_OMEGA:
        continue
      # Score tie-breakers:
      # prefer candidate that causes larger assigned color,
      # then prefer smaller new_omega, then shorter interval (less risk),
      # then lexicographically earlier placement for determinism.
      score = (assigned, -new_omega, - (cand[1] - cand[0]), -cand[0])
      if (assigned, -new_omega, - (cand[1] - cand[0])) > (best_score[0], best_score[1], 0):
        best_score = (assigned, -new_omega, cand)

    # If no candidate respects omega, fallback to inserting a short filler that doesn't change omega
    chosen = None
    if best_score[2] is not None:
      chosen = best_score[2]
    else:
      # create a tiny disjoint interval far to the right to avoid overlap growth
      # place it after current max endpoint + 2
      if intervals:
        max_r = max(r for _, r in intervals)
      else:
        max_r = 0
      chosen = (max_r + 2, max_r + 3)

    # Append the chosen interval
    intervals.append(chosen)

    # Occasionally (every 50 steps) do a small compaction:
    # move intervals to positive range (translate) to keep numbers tidy
    if (step + 1) % 50 == 0:
      # Rebase so smallest left is 0
      lo = min(l for l, r in intervals)
      if lo != 0:
        delta = -lo
        intervals = [(l + delta, r + delta) for (l, r) in intervals]

    # Stop early if we have reached a fairly high FirstFit color pressure
    colors = _firstfit_colors(intervals)
    if max(colors) >= 40 and len(intervals) >= 400:
      # safety stop: we already forced many colors with bounded omega
      break

  # Final ensure length within CAP
  if len(intervals) > CAP:
    intervals = intervals[:CAP]

  # Final normalization: translate so minimal left endpoint is 0 (clean output)
  if intervals:
    lo = min(l for l, r in intervals)
    if lo < 0:
      delta = -lo
      intervals = [(l + delta, r + delta) for (l, r) in intervals]

  return intervals

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()