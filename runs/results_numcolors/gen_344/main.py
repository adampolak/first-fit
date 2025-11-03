# EVOLVE-BLOCK-START

def construct_intervals(enable_alt_microphase=True):
  """
  Fractal-wave, region-fragmented interval generator.
  A deterministic, non-KT-builder designed to stress FirstFit via
  region-local density patterns rather than a multi-template spine.

  Returns:
    intervals: list of (l, r) pairs representing open intervals.
  """
  CAP = 9800  # maximum allowed intervals (as in prior runs)
  # Axis setup
  BLOCKS = 150           # number of non-overlapping blocks along the axis
  BLOCK_SPAN = 1_000_000.0 # span per block (units)
  MIN_WIDTH = 16.0
  MAX_WIDTH = 4000.0

  # Tiny deterministic PRNG (xorshift-ish) to seed local randomness
  seed = 0x9E3779B9
  def prng():
    nonlocal seed
    while True:
      x = seed & 0xFFFFFFFF
      x ^= (x << 13) & 0xFFFFFFFF
      x ^= (x >> 17) & 0xFFFFFFFF
      x ^= (x << 5) & 0xFFFFFFFF
      seed = x
      yield x

  RNG = prng()

  # Helper: emit a micro-interval inside a block [L, R)
  def emit_in_block(L, R, t, w_hint):
    # center inside block
    c = (L + R) * 0.5
    # width bounded, with deterministic jitter
    j = (next(RNG) & 0xF) / 16.0
    w = max(MIN_WIDTH, min(MAX_WIDTH, w_hint * (0.8 + j)))
    half = w * 0.5
    l = max(L, c - half)
    r = min(R, c + half)
    if r <= l:
      r = l + 1e-3
    return (l, r)

  # Build axis blocks
  blocks = []
  for b in range(BLOCKS):
    L = b * BLOCK_SPAN
    R = L + BLOCK_SPAN
    blocks.append((L, R))

  intervals = []

  # Stage 1: per-block micro-density
  # Emit 1-3 intervals per block, overlapping the same backbone region (the block's span)
  # to create strong local pressure while keeping omega modest by not over-densifying globally.
  for b, (L, R) in enumerate(blocks):
    if len(intervals) >= CAP:
      break
    m = 1 + ((b * 7) % 3)  # 1,2,3 in a deterministic rotating fashion
    for t in range(m):
      if len(intervals) >= CAP:
        break
      # width hint varies per block to diversify density
      w_hint =  MIN_WIDTH + ((b * 13 + t * 37) %  (MAX_WIDTH - MIN_WIDTH))
      intervals.append(emit_in_block(L, R, t, w_hint))

  # Stage 2: cross-block micro-extensions
  # A second pass across blocks, adding a small second interval per block but
  # offsetting them so that not all cross the exact same center. Keeps omega controlled.
  for b, (L, R) in enumerate(blocks):
    if len(intervals) >= CAP:
      break
    c = (L + R) * 0.5
    w = MIN_WIDTH * 0.75
    l = max(L, c - w)
    r = min(R, c + w)
    if r > l:
      intervals.append((l, r))

  # Stage 3: a handful of longer connectors spanning multiple blocks
  # These are appended to introduce cross-scale interactions without exploding omega.
  # They are chosen deterministically to ensure reproducibility.
  connectors_to_add = 60
  centers = [ (L+R)/2.0 for (L,R) in blocks ]
  for i in range(connectors_to_add):
    if len(intervals) >= CAP:
      break
    idx = (i * 11) % len(centers)
    center = centers[idx]
    span_blocks = max(1, BLOCKS // 4)
    span_len = span_blocks * BLOCK_SPAN / 4.0
    l = center - span_len
    r = l + max(100.0, span_len * 0.75)
    if r > l:
      intervals.append((l, r))

  # Normalize and finalize
  if not intervals:
    return []

  # Ensure endpoints are non-negative
  min_l = min(l for l, _ in intervals)
  if min_l < 0:
    intervals = [(l - min_l, r - min_l) for (l, r) in intervals]

  # Cap to CAP intervals
  if len(intervals) > CAP:
    intervals = intervals[:CAP]

  # Round endpoints for reproducibility
  final = []
  for (l, r) in intervals:
    if r <= l:
      r = l + 1e-3
    final.append((round(l, 6), round(r, 6)))
  return final

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()