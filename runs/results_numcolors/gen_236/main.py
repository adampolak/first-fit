# EVOLVE-BLOCK-START

def construct_intervals(rounds=6,
                        rotate_starts=True,
                        reverse_block_parity=True,
                        interleave_blocks=True,
                        phase2_iters=1,
                        cross4_enabled=True):
  """
  Rotating-template KT spine with dual micro-phases and long-range connectors.
  Returns intervals (l, r) in FirstFit presentation order.
  """
  CAP = 9800
  CAP_MARGIN = 32
  BASE_SEED = 0xC0FFEE123456789

  # Four strong KT templates
  template_bank = [
    (2, 6, 10, 14),
    (1, 5, 9, 13),
    (3, 7, 11, 15),
    (4, 8, 12, 16),
  ]

  # Seed with a single unit interval
  T = [(0, 1)]

  # Span helper
  def _span(T):
    lo = min(l for l, r in T)
    hi = max(r for l, r in T)
    d = hi - lo
    return lo, hi, d if d > 0 else 1

  # KT-round assembler
  def _apply_round(T_cur, starts, do_inter, rev_order, add_x4):
    lo, hi, d = _span(T_cur)
    blocks = []
    for s in starts:
      base = s * d - lo
      blocks.append([(l+base, r+base) for (l, r) in T_cur])
    S = []
    if do_inter:
      idxs = list(range(4))
      if rev_order: idxs.reverse()
      ml = max(len(b) for b in blocks)
      for i in range(ml):
        for j in idxs:
          if i < len(blocks[j]):
            S.append(blocks[j][i])
    else:
      if rev_order:
        blocks = list(reversed(blocks))
      for b in blocks:
        S.extend(b)
    # four KT connectors
    s0, s1, s2, s3 = starts
    S.append(((s0-1)*d, (s1-1)*d))
    S.append(((s2+2)*d, (s3+2)*d))
    S.append(((s0+2)*d, (s2-1)*d))
    S.append(((s1+2)*d, (s3-1)*d))
    # optional cross4
    if add_x4:
      S.append(((s0+4)*d, (s3+4)*d))
    return S

  # Stage 1: KT spine
  for ridx in range(rounds):
    # estimate next size
    nxt_conn = 4 + (1 if cross4_enabled else 0)
    if 4*len(T) + nxt_conn > CAP:
      break
    # deterministic seed per round
    seed = (BASE_SEED ^ (ridx * 0x9E3779B97F4A7C15)) & ((1<<64)-1)
    # choose starts
    if rotate_starts:
      starts = template_bank[seed % len(template_bank)]
    else:
      starts = template_bank[0]
    do_inter = interleave_blocks and ((seed>>1)&1==0)
    rev_ord = reverse_block_parity and ((seed>>2)&1==1)
    add_x4 = cross4_enabled and ((seed>>3)&1==1)
    T = _apply_round(T, starts, do_inter, rev_ord, add_x4)
    if len(T) >= CAP:
      return T[:CAP]
  if len(T) >= CAP - CAP_MARGIN:
    return T[:CAP]

  lo, hi, d = _span(T)
  # Micro-phase A: three tail caps
  def cap_at(a,b):
    L = lo + max(1,int(round(a*d)))
    R = lo + max(1,int(round(b*d)))
    return (L, R if R>L else L+1)
  for a,b in [(0.05,0.60),(0.25,0.75),(0.75,0.92)]:
    if len(T)>=CAP: break
    T.append(cap_at(a,b))
  if len(T) >= CAP - CAP_MARGIN:
    return T[:CAP]

  # micro-builder
  def build_micro(T_cur, rid, budget, windows):
    if budget <= 8 or not T_cur:
      return []
    glo, ghi, G = _span(T_cur)
    # thin seed
    n = len(T_cur)
    seed_sz = max(8, min(40, n//260))
    stride = max(1, n//seed_sz)
    off = (BASE_SEED ^ rid) % stride
    U = [T_cur[(i+off)%n] for i in range(0,n,stride)][:seed_sz]
    if not U: return []
    ulo = min(l for l,_ in U)
    # build blocks
    blocks = []
    for fa, fb in windows:
      win = glo + int(round(fa*G))
      base = win - ulo
      blk = [(l+base,r+base) for (l,r) in U]
      if int((fa+rid)*100)%2==0:
        blk.reverse()
      blocks.append(blk)
    # interleave
    micro = []
    ml = max(len(b) for b in blocks)
    order = list(range(len(blocks)))
    if rid&1: order.reverse()
    for i in range(ml):
      for j in order:
        if i < len(blocks[j]):
          micro.append(blocks[j][i])
    # connectors
    for fa, fb in [(0.08,0.30),(0.60,0.92),(0.26,0.56),(0.44,0.78)]:
      a = glo+int(round(fa*G))
      b = glo+int(round(fb*G))
      if b>a:
        micro.append((a,b))
    return micro[:budget]

  # Stage 2: dual micro-phases
  windows1 = [(0.12,0.22),(0.35,0.45),(0.58,0.68),(0.80,0.90)]
  windows2 = [(0.05,0.15),(0.28,0.38),(0.60,0.70),(0.82,0.92)]
  for rid in range(min(phase2_iters,2)):
    if len(T) > CAP - CAP_MARGIN: break
    room = CAP - len(T)
    win = windows1 if rid==0 else windows2
    mic = build_micro(T, rid, room, win)
    if not mic: break
    T.extend(mic)
  if len(T) >= CAP - CAP_MARGIN:
    return T[:CAP]

  # Final long-range cross-scale connectors
  lo, hi, d = _span(T)
  cross = [
    (lo+int(0.02*d), hi-int(0.02*d)),
    (lo+int(0.15*d), hi-int(0.15*d)),
    (lo+int(0.30*d), hi-int(0.30*d)),
    (lo+int(0.45*d), hi-int(0.45*d)),
    (lo+int(0.55*d), hi-int(0.55*d)),
  ]
  for c in cross:
    if len(T)>=CAP: break
    if c[1]>c[0]:
      T.append(c)

  return T[:CAP]

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()