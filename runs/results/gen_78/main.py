# EVOLVE-BLOCK-START

import itertools

# ---------- Evaluator ----------
def firstfit_colors(intervals):
    last_end = []
    for l, r in intervals:
        placed = False
        for i, le in enumerate(last_end):
            if l >= le:
                last_end[i] = r
                placed = True
                break
        if not placed:
            last_end.append(r)
    return len(last_end)

def clique_number(intervals):
    events = []
    for l, r in intervals:
        if l < r:
            events.append((l, +1))
            events.append((r, -1))
    events.sort(key=lambda e: (e[0], e[1]))
    cur = best = 0
    for _, delta in events:
        cur += delta
        if cur > best:
            best = cur
    return best

def normalize_intervals(intervals):
    if not intervals:
        return []
    pts = sorted({x for seg in intervals for x in seg})
    mapping = {x:i*2 for i,x in enumerate(pts)}
    return [(mapping[l], mapping[r]) for l, r in intervals]

# ---------- Pattern Generator ----------
class PatternGen:
    @staticmethod
    def generate(base_seed, depth, offsets, blockers, extra_first):
        T = list(base_seed)
        for lvl in range(depth):
            lo = min(l for l,_ in T)
            hi = max(r for _,r in T)
            d = hi - lo
            offs = list(offsets)
            if extra_first and lvl == 0:
                offs.append(offsets[-1] + 4)
            S = []
            for o in offs:
                off = d * o - lo
                S.extend([(l+off, r+off) for l,r in T])
            for a, b in blockers:
                S.append((d*a, d*b))
            T = S
        return T

# ---------- Search Engine ----------
class SearchEngine:
    def __init__(self):
        self.depths = [3,4,5]
        self.offset_sets = [
            (2,6,10,14),
            (1,5,9,13),
            (3,7,11,15),
            (0,4,8,12),
        ]
        self.blocker_sets = [
            [(1,5),(12,16),(4,9),(8,13)],
            [(1,6),(11,16),(3,9),(7,13)],
            [(2,6),(12,16),(4,9),(8,13)],
        ]
        self.extra_choices = [False, True]
        self.base_seed = [(0.0, 1.0)]

    def find_best(self):
        best = None  # (ratio, cols, om, raw_intervals)
        for depth, offs, blks, extra in itertools.product(
            self.depths, self.offset_sets, self.blocker_sets, self.extra_choices
        ):
            raw = PatternGen.generate(self.base_seed, depth, offs, blks, extra)
            norm = normalize_intervals(raw)
            om = clique_number(norm)
            if om == 0:
                continue
            cols = firstfit_colors(norm)
            ratio = cols/om
            if best is None or ratio > best[0] + 1e-12 \
               or (abs(ratio-best[0])<=1e-12 and len(norm) < len(best[3])):
                best = (ratio, cols, om, raw)
        if best is None:
            # fallback to canonical depth=4
            raw = PatternGen.generate(self.base_seed, 4, self.offset_sets[0], self.blocker_sets[0], False)
            return normalize_intervals(raw), clique_number(normalize_intervals(raw))
        return normalize_intervals(best[3]), best[2]

# ---------- Augmentor ----------
class Augmentor:
    @staticmethod
    def inject_waves(intervals, target_omega, wave_lengths=(2,4), max_waves=60):
        cur = list(intervals)
        base_cols = firstfit_colors(cur)
        waves_added = 0
        for L in wave_lengths:
            if waves_added >= max_waves:
                break
            min_x = min(l for l,_ in cur)
            max_x = max(r for _,r in cur)
            for x in range(min_x, max_x - L + 1):
                if waves_added >= max_waves:
                    break
                candidate = (x, x+L)
                if clique_number(cur + [candidate]) > target_omega:
                    continue
                new_cols = firstfit_colors(cur + [candidate])
                if new_cols > base_cols:
                    cur.append(candidate)
                    base_cols = new_cols
                    waves_added += 1
        return cur

# ---------- Pruner ----------
class Pruner:
    @staticmethod
    def prune(intervals, keep_ratio):
        T = list(intervals)
        changed = True
        while changed:
            changed = False
            # remove in descending length order
            order = sorted(range(len(T)), key=lambda i: -(T[i][1]-T[i][0]))
            for idx in order:
                cand = T[:idx] + T[idx+1:]
                norm = normalize_intervals(cand)
                om = clique_number(norm)
                if om == 0:
                    continue
                cols = firstfit_colors(norm)
                if cols/om >= keep_ratio - 1e-12:
                    T = cand
                    changed = True
                    break
        return T

# ---------- Main Construction ----------
def construct_intervals():
    # 1) Search over parameter space
    engine = SearchEngine()
    base_norm, base_om = engine.find_best()

    # 2) Augment with waves to push FirstFit
    augmented = Augmentor.inject_waves(base_norm, base_om)

    # 3) Prune redundancies while preserving ratio
    final_ratio = firstfit_colors(augmented) / max(1, clique_number(augmented))
    pruned = Pruner.prune(augmented, final_ratio)

    # 4) Final normalization
    return normalize_intervals(pruned)

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()