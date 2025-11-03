# EVOLVE-BLOCK-START

from typing import List, Tuple

Interval = Tuple[int, int]


class SeedChain:
  """Deterministic integer generator; no true randomness to preserve reproducibility."""
  def __init__(self, base: int = 1337):
    self.state = (base ^ 0x9e3779b97f4a7c15) & ((1 << 64) - 1)

  def next(self) -> int:
    # XorShift-like step, deterministic
    x = self.state
    x ^= (x << 13) & ((1 << 64) - 1)
    x ^= (x >> 7)
    x ^= (x << 17) & ((1 << 64) - 1)
    self.state = x & ((1 << 64) - 1)
    return self.state

  def pick(self, lo: int, hi: int) -> int:
    if hi <= lo:
      return lo
    return lo + (self.next() % (hi - lo + 1))


class CapManager:
  """Tracks capacity and vetoes additions that would exceed CAP."""
  def __init__(self, cap: int):
    self.CAP = cap

  def room(self, n: int) -> int:
    return max(0, self.CAP - n)

  def can_add(self, n: int, add: int) -> bool:
    return n + add <= self.CAP

  def trim(self, seq: List[Interval], extra: List[Interval]) -> List[Interval]:
    room = self.room(len(seq))
    return extra[:room]


def _span_delta(T: List[Interval]) -> Tuple[int, int, int]:
  lo = min(l for l, _ in T)
  hi = max(r for _, r in T)
  delta = hi - lo
  if delta <= 0:
    delta = 1
  return lo, hi, delta


def _translate_block(T: List[Interval], base_shift: int) -> List[Interval]:
  return [(l + base_shift, r + base_shift) for (l, r) in T]


def _assemble_blocks(blocks: List[List[Interval]], interleave: bool, reverse_order: bool, rot: int = 0) -> List[Interval]:
  if not blocks:
    return []
  k = len(blocks)
  order = list(range(k))
  if reverse_order:
    order.reverse()
  if k > 0 and rot % k != 0:
    r = rot % k
    order = order[r:] + order[:r]

  if not interleave:
    S = []
    for idx in order:
      S.extend(blocks[idx])
    return S

  S = []
  maxlen = max(len(b) for b in blocks)
  for i in range(maxlen):
    for idx in order:
      blk = blocks[idx]
      if i < len(blk):
        S.append(blk[i])
  return S


def _classic_connectors(starts: Tuple[int, int, int, int], delta: int) -> List[Interval]:
  s0, s1, s2, s3 = starts
  return [
    ((s0 - 1) * delta, (s1 - 1) * delta),
    ((s2 + 2) * delta, (s3 + 2) * delta),
    ((s0 + 2) * delta, (s2 - 1) * delta),
    ((s1 + 2) * delta, (s3 - 1) * delta),
  ]


def _long_range_connectors(starts: Tuple[int, int, int, int], delta: int, mode: int = 0) -> List[Interval]:
  # Conservative long-range connectors spanning approx 4–6 deltas.
  # Anchored to spine starts; designed to couple colors without inflating omega drastically.
  s0, s1, s2, s3 = starts
  L = []
  if mode % 2 == 0:
    # Slightly inside edges to avoid peak density points
    L.append(((s0 + 1) * delta, (s2 + 5) * delta))
    L.append(((s1 + 2) * delta, (s3 + 6) * delta))
  else:
    # Crossed pairing
    L.append(((s0 - 2) * delta, (s3 + 3) * delta))
    L.append(((s1 - 1) * delta, (s2 + 4) * delta))
  # Ensure valid
  out = []
  for a, b in L:
    if b > a:
      out.append((a, b))
  return out


def _thin_seed(current_T: List[Interval], max_seed: int) -> List[Interval]:
  n = len(current_T)
  if n == 0 or max_seed <= 0:
    return []
  step = max(1, n // max_seed)
  return current_T[::step][:max_seed]


def _build_micro_windows(current_T: List[Interval],
                         seed_sz: int,
                         window_fracs: List[Tuple[float, float]],
                         interleave: bool,
                         reverse_order: bool,
                         rot: int,
                         add_connectors: bool) -> List[Interval]:
  if not current_T:
    return []
  glo = min(l for l, _ in current_T)
  ghi = max(r for _, r in current_T)
  G = max(1, ghi - glo)

  U = _thin_seed(current_T, seed_sz)
  if not U:
    return []
  ulo = min(l for l, _ in U)

  blocks = []
  for (fa, fb) in window_fracs:
    fa = max(0.05, min(0.90, fa))
    fb = max(0.10, min(0.95, fb))
    win_lo = glo + int(round(fa * G))
    base = win_lo - ulo
    block = _translate_block(U, base)
    # Weak internal symmetry break
    if ((int(round(fa * 100)) // 5) + rot) % 2 == 1:
      block = list(reversed(block))
    blocks.append(block)

  micro = _assemble_blocks(blocks, interleave=interleave, reverse_order=reverse_order, rot=rot)

  if add_connectors:
    connectors = [
      (glo + int(round(0.08 * G)), glo + int(round(0.30 * G))),
      (glo + int(round(0.26 * G)), glo + int(round(0.56 * G))),
      (glo + int(round(0.44 * G)), glo + int(round(0.78 * G))),
      (glo + int(round(0.60 * G)), glo + int(round(0.92 * G))),
    ]
    for a, b in connectors:
      if b > a:
        micro.append((a, b))

  return micro


def construct_intervals(seed_count: int = 1) -> List[Interval]:
  """
  Construct a sequence of intervals on the real line, in the order presented to FirstFit,
  engineered to drive up FirstFit's color usage while keeping omega modest (<= ~10).

  Input:
    seed_count: optional integer steering the deterministic seed chain.
  Output:
    List of (l, r) integer open intervals with r > l, total count <= 9800.
  """
  CAP = 9800
  capman = CapManager(CAP)
  BASE_SEED = int(seed_count) if isinstance(seed_count, int) else 1
  rng = SeedChain(base=0xC0FFEE ^ (BASE_SEED & 0xffffffff))

  # Spine templates (rotating KT-style starts)
  template_bank = [
    (2, 6, 10, 14),
    (1, 5, 9, 13),
    (3, 7, 11, 15),
    (4, 8, 12, 16),
  ]

  # Seed: one unit interval (multi-seed path retained but conservative)
  if seed_count and seed_count > 1:
    seeds = min(4, int(seed_count))
    step = 3
    T: List[Interval] = [(i * step, i * step + 1) for i in range(seeds)]
  else:
    T: List[Interval] = [(0, 1)]

  # Helper: apply one spine round with optional CAP-aware densification K
  def apply_spine_round(current_T: List[Interval], starts: Tuple[int, int, int, int],
                        round_idx: int, allow_densify: bool) -> List[Interval]:
    lo, hi, delta = _span_delta(current_T)

    # Decide assembly policy: even rounds interleave, odd rounds sequential reversed
    interleave = (round_idx % 2 == 0)
    reverse_order = (round_idx % 2 == 1)

    # Primary blocks (4 starts)
    blocks = []
    for s in starts:
      base = s * delta - lo
      blocks.append(_translate_block(current_T, base))

    # CAP-aware density multiplier K: inject one midpoint block between s1 and s2 on even rounds only
    # and only when capacity allows; deterministic based on seed chain and round index.
    if allow_densify and (round_idx % 2 == 0):
      s0, s1, s2, s3 = starts
      mid = (s1 + s2) // 2
      extra_base = mid * delta - lo
      extra_block = _translate_block(current_T, extra_base)
      # Only add if within CAP budget (estimate: len(current_T))
      if capman.can_add(len([]), 0) and capman.can_add(len(current_T) * 5, 0):  # trivial always True
        # Real guard: check if adding would exceed CAP after assembly (approximate)
        if capman.can_add(len(T), len(current_T)):
          blocks.insert(2, extra_block)  # place near middle to enhance coupling

    # Assemble blocks under parity schedule
    S = _assemble_blocks(blocks, interleave=interleave, reverse_order=reverse_order, rot=round_idx)

    # Classic connectors
    S.extend(_classic_connectors(starts, delta))

    # Deterministic long-range cross-scale connectors, sparsely added on last two rounds
    if round_idx >= 4:
      mode = rng.pick(0, 1)
      S.extend(_long_range_connectors(starts, delta, mode=mode))

    return S

  # Stage 1: backbone rounds (cap-checked), targeting up to six rounds
  max_rounds = 6
  for ridx in range(max_rounds):
    # Predict next size: 4*size + 4 (rough baseline without densify)
    projected = 4 * len(T) + 4
    if projected > CAP:
      break
    starts = template_bank[ridx % len(template_bank)]
    allow_densify = (ridx in (2, 4))  # densify only on a couple of mid/late rounds
    T = apply_spine_round(T, starts, ridx, allow_densify=allow_densify)
    if len(T) >= CAP:
      T = T[:CAP]
      return T

  # Early return if near capacity
  if len(T) >= CAP - 8:
    return T[:CAP]

  # Micro-phase A: deterministic window family 1 with interleaving
  loA, hiA, deltaA = _span_delta(T)
  seed_sz_A = max(8, min(36, len(T) // 280))
  window_fracs_A = [(0.12, 0.22), (0.35, 0.45), (0.58, 0.68), (0.80, 0.90)]
  microA = _build_micro_windows(
    T, seed_sz=seed_sz_A, window_fracs=window_fracs_A,
    interleave=True, reverse_order=False, rot=BASE_SEED % 5, add_connectors=True
  )
  if microA:
    microA = capman.trim(T, microA)
    T.extend(microA)

  if len(T) >= CAP - 8:
    return T[:CAP]

  # Micro-phase B: deterministic window family 2 with reverse order
  loB, hiB, deltaB = _span_delta(T)
  seed_sz_B = max(8, min(40, len(T) // 250))
  window_fracs_B = [(0.06, 0.14), (0.28, 0.38), (0.54, 0.64), (0.78, 0.88)]
  microB = _build_micro_windows(
    T, seed_sz=seed_sz_B, window_fracs=window_fracs_B,
    interleave=True, reverse_order=True, rot=(BASE_SEED + 1) % 7, add_connectors=True
  )
  if microB:
    microB = capman.trim(T, microB)
    T.extend(microB)

  # Final CAP guard
  if len(T) > CAP:
    T = T[:CAP]

  # Normalize: ensure non-negative integer endpoints and r > l
  if not T:
    return []
  min_l = min(l for l, _ in T)
  if min_l < 0:
    T = [(l - min_l, r - min_l) for (l, r) in T]

  out: List[Interval] = []
  for (l, r) in T:
    li = int(l)
    ri = int(r)
    if ri <= li:
      ri = li + 1
    out.append((li, ri))

  if len(out) > CAP:
    out = out[:CAP]
  return out


def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals(**kwargs)

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()