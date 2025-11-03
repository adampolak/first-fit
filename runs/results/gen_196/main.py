# EVOLVE-BLOCK-START
from functools import lru_cache
from math import gcd
from itertools import product

# ------------------------------
# Evaluator module
# ------------------------------
class Evaluator:
    @staticmethod
    def normalize(intervals):
        """Map endpoints to a compact even-integer grid."""
        if not intervals:
            return []
        pts = sorted({x for l, r in intervals for x in (l, r)})
        m = {p: 2*i for i, p in enumerate(pts)}
        return [(m[l], m[r]) for l, r in intervals]

    @staticmethod
    def clique_number(intervals):
        """Sweep-line for open intervals."""
        ev = []
        for l, r in intervals:
            if l < r:
                ev.append((l, +1))
                ev.append((r, -1))
        if not ev:
            return 0
        ev.sort(key=lambda e: (e[0], 0 if e[1] == -1 else 1))
        cur = best = 0
        for _, d in ev:
            cur += d
            if cur > best:
                best = cur
        return best

    @staticmethod
    def firstfit_colors(intervals):
        """FirstFit simulation via per-color last_end tracking."""
        last_end = []
        for l, r in intervals:
            placed = False
            for i, e in enumerate(last_end):
                if l >= e:
                    last_end[i] = r
                    placed = True
                    break
            if not placed:
                last_end.append(r)
        return len(last_end)

    @lru_cache(maxsize=None)
    def evaluate_raw(self, raw_tuple):
        """
        Evaluate a raw pattern (tuple of tuples). Returns
        (score, omega, colors, n, normalized_intervals).
        """
        raw = [tuple(iv) for iv in raw_tuple]
        Tn = self.normalize(raw)
        n = len(Tn)
        if n == 0:
            return (-1.0, 0, 0, n, Tn)
        om = self.clique_number(Tn)
        if om == 0:
            return (-1.0, 0, 0, n, Tn)
        cols = self.firstfit_colors(Tn)
        ratio = cols / om
        # small penalty for size
        score = ratio - 1e-6 * (n/1000.0)
        return (score, om, cols, n, Tn)

# ------------------------------
# Pruner module
# ------------------------------
class Pruner:
    @staticmethod
    def shrink_strict(raw, target_cols, target_om, max_passes=3):
        """Greedy removal keeping exact (cols,om)."""
        ev = Evaluator()
        cur = list(raw)
        for _ in range(max_passes):
            removed = False
            # try removing largest intervals first
            order = sorted(range(len(cur)),
                           key=lambda i: -(cur[i][1]-cur[i][0]))
            for i in order:
                cand = tuple(cur[:i] + cur[i+1:])
                score, om, cols, n, _ = ev.evaluate_raw(cand)
                if om == target_om and cols == target_cols:
                    cur = list(cand)
                    removed = True
                    break
            if not removed:
                break
        return tuple(cur)

    @staticmethod
    def final_shrink(raw, target_cols, target_om):
        """Single-interval removals preserving ratio and ω."""
        ev = Evaluator()
        cur = list(raw)
        while True:
            removed = False
            for i in range(len(cur)):
                cand = tuple(cur[:i] + cur[i+1:])
                score, om, cols, n, _ = ev.evaluate_raw(cand)
                if om == target_om and cols == target_cols:
                    cur = list(cand)
                    removed = True
                    break
            if not removed:
                break
        return tuple(cur)

# ------------------------------
# Pattern builder module
# ------------------------------
class PatternBuilder:
    @staticmethod
    @lru_cache(maxsize=None)
    def build_raw(depth, offsets, blockers, translation, anchor, wave, perturb_idx, perturb_delta):
        """
        Build a raw tuple-of-tuples pattern with:
          depth: recursion depth
          offsets: 4-tuple of floats
          blockers: tuple of (a,b) floats
          translation: 'left' or 'center'
          anchor: 'left' or 'center'
          wave: 'none','short','long'
          perturb_idx: which offset to perturb (0..3) or -1 for none
          perturb_delta: small float
        """
        T = [(0.0, 1.0)]
        for lvl in range(depth):
            lo = min(l for l, r in T)
            hi = max(r for l, r in T)
            delta = hi - lo
            center = (lo+hi)/2.0
            # apply perturb
            offs = []
            for idx, s in enumerate(offsets):
                s2 = s + (perturb_delta if idx==perturb_idx else 0.0)
                offs.append(s2)
            # make copies
            S = []
            for s in offs:
                off = (delta*s - (lo if translation=='left' else center))
                for l,r in T:
                    S.append((l+off, r+off))
            # add blockers
            for a,b in blockers:
                offa = delta*a - (0.0 if anchor=='left' else center)
                offb = delta*b - (0.0 if anchor=='left' else center)
                S.append((offa, offb))
            T = S
        # wave injection
        if wave!='none' and T:
            lo = min(l for l,r in T)
            hi = max(r for l,r in T)
            delta = hi-lo
            if wave=='short':
                w = delta*0.05
                pos = lo + delta*1.5
                T.append((pos, pos+w))
            else:  # 'long'
                w = delta*0.2
                pos = lo + delta*3.0
                T.append((pos, pos+w))
        # return as hashable tuple
        return tuple(T)

# ------------------------------
# Main construction
# ------------------------------
def construct_intervals():
    ev = Evaluator()
    best = None  # (score, raw_tuple, om, cols)
    # search parameters
    offsets_list = [
        (2.0,6.0,10.0,14.0),
        (1.0,5.0,9.0,13.0),
        (3.0,7.0,11.0,15.0),
        (0.0,4.0,8.0,12.0),
    ]
    blockers_list = [
        ((1.0,5.0),(12.0,16.0),(4.0,9.0),(8.0,13.0)),
        ((0.0,4.0),(6.0,10.0),(8.0,12.0),(14.0,18.0)),
        ((2.0,6.0),(4.0,8.0),(10.0,14.0),(12.0,16.0)),
    ]
    translations = ['left','center']
    anchors = ['left','center']
    depths = [3,4,5]
    waves = ['none','short','long']
    perturbs = [-1,0,1]  # -1 means no perturb, 0..3 index into offsets
    perturb_delta = 0.1

    # enumerate
    for (depth, offsets, blockers,
         translation, anchor, wave, pidx) in product(
            depths, offsets_list, blockers_list,
            translations, anchors, waves, perturbs):
        raw = PatternBuilder.build_raw(
            depth, offsets, blockers,
            translation, anchor,
            wave, pidx-1, perturb_delta
        )
        score, om, cols, n, _ = ev.evaluate_raw(raw)
        if om==0:
            continue
        if best is None or score > best[0] + 1e-9:
            best = (score, raw, om, cols)
        elif best and abs(score-best[0])<=1e-9:
            # tie-break fewer intervals then more colors
            _, _, bom, bcols = best
            curr_n = len(raw)
            if curr_n < len(best[1]) or (curr_n==len(best[1]) and cols>bcols):
                best = (score, raw, om, cols)

    # fallback
    if best is None:
        best_raw = PatternBuilder.build_raw(4, offsets_list[0], blockers_list[0],
                                           'left','left','none',-1,0.0)
        score, om, cols, n, norm = ev.evaluate_raw(best_raw)
        return norm

    _, best_raw, best_om, best_cols = best
    # prune
    strict = Pruner.shrink_strict(best_raw, best_cols, best_om, max_passes=4)
    final = Pruner.final_shrink(strict, best_cols, best_om)
    # normalized
    _, om2, cols2, n2, norm = ev.evaluate_raw(final)
    # sanity
    if om2==best_om and cols2==best_cols:
        return norm
    # else fallback
    _, _, _, _, norm0 = ev.evaluate_raw(best_raw)
    return norm0

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()