# EVOLVE-BLOCK-START

def _span_delta(T):
  lo = min(l for l, _ in T)
  hi = max(r for _, r in T)
  delta = hi - lo
  if delta <= 0:
    delta = 1
  return lo, hi, delta

def _build_blocks_sequential(T, starts, delta, lo, K=1):
  """
  Build four translated blocks in a simple, sequential (non-interleaved) fashion.
  This ordering is known to push FirstFit colors higher while keeping omega in check.
  """
  S = []
  for s in starts:
    base = s * delta * K - lo
    for (l, r) in T:
      S.append((l + base,   r + base))
  return S

def _append_connectors(S, starts, delta, add_cross4=False):
  """
  Append the classic four connectors. Optionally add a single long-range cross4
  (final-round only) to enhance pressure with minimal omega impact.
  """
  s0, s1, s2, s3 = starts
  S.append(((s0 - 1) * delta, (s1 - 1) * delta))  # left cap
  S.append(((s2 + 2) * delta, (s3 + 2) * delta))  # right cap
  S.append(((s0 + 2) * delta, (s2 - 1) * delta))  # cross 1
  S.append(((s1 + 2) * delta, (s3 - 1) * delta))  # cross 2
  if add_cross4:
    S.append(((s0 + 4) * delta, (s3 + 4) * delta))  # long-range cross (final-round)

def _append_micro_tail(T, max_extra=8):
  """
  Append a tiny tail of long caps near the end of T to boost FF late with modest omega impact.
  Caps are inserted near the tail to intersect many active colors.
  """
  if not T or max_extra <= 0:
    return T
  lo, hi, delta = _span_delta(T)
  d2 = max(1, delta // 4)
  micro = [
      (lo + 1 * d2, lo + 5 * d2),
      (lo + 3 * d2, lo + 8 * d2),
      (hi - 6 * d2, hi - 2 * d2),
      (hi - 8 * d2, hi - 3 * d2),
  ]
  for i, interval in enumerate(micro[:max_extra]):
    pos = len(T) - (i * 2 + 1)
    if pos < 0:
      T.append(interval)
    else:
      T.insert(pos, interval)
  return T

def construct_intervals(seed_count=1):
  """
  Construct a deterministic interval sequence with a KT spine and a small tail,
  designed to maximize FirstFit colors while keeping clique size (omega) small.
  """
  CAP = 9800

  # Four-template spine bank (rotated round-robin)
  template_bank = [
    (2, 6, 10, 14),  # T1 classic KT
    (1, 5, 9, 13),   # T2 left-shifted
    (3, 7, 11, 15),  # T3 right-shifted
    (4, 8, 12, 16),  # T4 stretched-right
  ]

  # Seed: single unit interval
  T = [(0, 1)]

  rounds = 6
  for ridx in range(rounds):
    lo, hi, delta = _span_delta(T)
    starts = template_bank[ridx % len(template_bank)]
    S = _build_blocks_sequential(T, starts, delta, lo, K=1)
    _append_connectors(S, starts, delta, add_cross4=(ridx == rounds - 1))
    T = S
    if len(T) >= CAP:
      T = T[:CAP]
      return T

  # Tiny micro tail to boost late FF interactions, cap-bounded
  room = CAP - len(T)
  if room > 0:
    T = _append_micro_tail(T, max_extra=min(8, room))

  return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()