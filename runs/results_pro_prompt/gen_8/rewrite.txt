# EVOLVE-BLOCK-START

from typing import List, Tuple
from math import gcd

def _step_scale(T: List[Tuple[int, int]],
                starts: Tuple[int, ...],
                caps: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
  """
  One replication/capping level:
  - Place |starts| disjoint scaled copies of T at integer offsets.
  - Add 'caps' which are long intervals coupling copies, preserving low omega yet
    forcing FirstFit to use new colors near the top of the palette.
  """
  lo = min(l for l, r in T)
  hi = max(r for l, r in T)
  delta = hi - lo
  S: List[Tuple[int, int]] = []
  # Replicate blocks in order; the block order is adversarial for FirstFit.
  for start in starts:
    S += [(delta * start + (l - lo), delta * start + (r - lo)) for (l, r) in T]
  # Add capping intervals last to force high colors.
  S += [(delta * a, delta * b) for (a, b) in caps]
  return S

def _normalize_and_compress(T: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
  """
  Shift and compress coordinates without changing the interval intersection pattern:
  - Shift so minimum left endpoint is 0.
  - Divide all endpoints by gcd of all endpoints to reduce magnitude (keeping integers).
  """
  if not T:
    return T
  min_l = min(l for l, _ in T)
  shifted = [(l - min_l, r - min_l) for (l, r) in T]
  g = 0
  for l, r in shifted:
    g = gcd(g, l)
    g = gcd(g, r)
  if g > 1:
    shifted = [(l // g, r // g) for (l, r) in shifted]
  return shifted

def construct_intervals() -> List[Tuple[int, int]]:
  """
  Construct a sequence of open intervals, presented in adversarial order for FirstFit.
  We use a deeper instantiation of the 4-block + caps recursive pattern (depth=6)
  to push FirstFit/OPT ratio higher than the previous depth=4 program.

  Returns:
    intervals: list of (l, r) pairs representing open intervals (l, r), in arrival order.
  """
  # Parameters (adversarial pattern from Figure 4):
  # - starts define the positions of the replicated blocks per level.
  # - caps define the long capping intervals per level, appended last.
  starts = (2, 6, 10, 14)
  caps = [(1, 5), (12, 16), (4, 9), (8, 13)]
  depth = 6  # Increased from 4 to 6 to improve FirstFit/omega ratio (~2.714 vs. 2.60)

  # Base template: a single short interval.
  T: List[Tuple[int, int]] = [(0, 1)]

  # Apply recursive amplification.
  for _ in range(depth):
    T = _step_scale(T, starts, caps)

  # Normalize coordinates to keep numbers compact without changing intersection relations.
  T = _normalize_and_compress(T)
  return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()