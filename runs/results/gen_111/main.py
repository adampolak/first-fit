# EVOLVE-BLOCK-START
from dataclasses import dataclass
from itertools import product
from math import gcd

@dataclass(frozen=True)
class PatternParams:
    depth: int
    offsets: tuple
    blockers: tuple
    translation: str   # 'left' or 'center'
    schedule: str      # 'after','before','split'

class PatternBuilder:
    @staticmethod
    def build(params: PatternParams):
        T = [(0.0, 1.0)]
        for _ in range(params.depth):
            lo = min(l for l, _ in T)
            hi = max(r for _, r in T)
            delta = hi - lo
            center = (lo + hi) / 2.0

            def make_copies():
                for start in params.offsets:
                    base = lo if params.translation == 'left' else center
                    offset = delta * start - base
                    for l, r in T:
                        yield (l + offset, r + offset)

            def make_blockers():
                for a, b in params.blockers:
                    if params.translation == 'left':
                        yield (delta * a, delta * b)
                    else:
                        yield (delta * a - center, delta * b - center)

            if params.schedule == 'before':
                T = list(make_blockers()) + list(make_copies())
            elif params.schedule == 'split':
                half = len(params.offsets) // 2
                first, second = params.offsets[:half], params.offsets[half:]
                def copies(starts):
                    base = lo if params.translation == 'left' else center
                    for start in starts:
                        offset = delta * start - base
                        for l, r in T:
                            yield (l + offset, r + offset)
                T = list(copies(first)) + list(make_blockers()) + list(copies(second))
            else:  # 'after'
                T = list(make_copies()) + list(make_blockers())
        return T

class Evaluator:
    @staticmethod
    def ff_count(intervals):
        last_end = []
        for l, r in intervals:
            for i, le in enumerate(last_end):
                if l >= le:
                    last_end[i] = r
                    break
            else:
                last_end.append(r)
        return len(last_end)

    @staticmethod
    def omega(intervals):
        events = []
        for l, r in intervals:
            if l < r:
                events.append((l, 1))
                events.append((r, -1))
        events.sort(key=lambda e: (e[0], 0 if e[1] < 0 else 1))
        cur = best = 0
        for _, t in events:
            cur += t
            if cur > best:
                best = cur
        return best

    @staticmethod
    def normalize(intervals):
        pts = sorted({x for seg in intervals for x in seg})
        coord = {x: 2*i for i, x in enumerate(pts)}
        norm = [(coord[l], coord[r]) for l, r in intervals]
        allv = [v for seg in norm for v in seg]
        g = 0
        for v in allv: 
            g = gcd(g, v)
        if g > 1:
            norm = [(l//g, r//g) for l, r in norm]
        min_l = min(l for l, _ in norm)
        if min_l < 0:
            norm = [(l - min_l, r - min_l) for l, r in norm]
        return norm

def construct_intervals():
    depths = [3, 4, 5]
    offset_sets = [(2,6,10,14), (1,5,9,13), (3,7,11,15), (0,4,8,12)]
    blocker_templates = [
        ((1,5),(12,16),(4,9),(8,13)),
        ((0,4),(6,10),(8,12),(14,18)),
        ((2,6),(4,8),(10,14),(12,16))
    ]
    translations = ['left','center']
    schedules = ['after','before','split']

    best = None  # (score_tuple, intervals)
    for depth, offs, blks, tr, sch in product(depths, offset_sets, blocker_templates, translations, schedules):
        params = PatternParams(depth, offs, blks, tr, sch)
        raw = PatternBuilder.build(params)
        n = len(raw)
        if n > 2000:
            continue
        om = Evaluator.omega(raw)
        if om <= 0:
            continue
        cols = Evaluator.ff_count(raw)
        ratio = cols / om
        score = (ratio, -n, cols)
        if best is None or score > best[0]:
            best = (score, raw)

    if best is None:
        # fallback to canonical depth=4
        params = PatternParams(4, (2,6,10,14), ((1,5),(12,16),(4,9),(8,13)), 'left', 'after')
        raw = PatternBuilder.build(params)
    else:
        raw = best[1]

    return Evaluator.normalize(raw)
# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()