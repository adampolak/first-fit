# EVOLVE-BLOCK-START

def construct_intervals():
    """
    Construct a sequence of open intervals that forces FirstFit to use many colors
    relative to the clique number (OPT). We use the classic 4-tiles + 4-blockers
    recursive pattern (cf. Fig. 4 in https://arxiv.org/abs/1506.00192), which
    grows the FirstFit colors roughly three times as fast as the clique, pushing
    the ratio toward 3 while keeping the instance moderately sized.

    Returns:
      intervals: list[(l, r)] of open intervals with integer coordinates.
    """
    # Configuration: keep instance size modest while improving the FF/OPT ratio.
    # With this cap, depth k=5 yields ~2388 intervals, FF=16, OPT=6 (≈2.6667).
    MAX_INTERVALS = 2500

    # Tiling positions and blocker windows (scaled each round by current span).
    tile_positions = (2, 6, 10, 14)
    blocker_windows = [(1, 5), (12, 16), (4, 9), (8, 13)]

    # Size growth: n_{k} = 4 n_{k-1} + 4, n_0 = 1  =>  n_k = ((7 * 4^k) - 4) / 3
    def size_for_depth(k: int) -> int:
        return ((7 * (4 ** k)) - 4) // 3

    # Choose the largest k such that size_for_depth(k) <= MAX_INTERVALS
    k = 0
    while size_for_depth(k + 1) <= MAX_INTERVALS:
        k += 1

    # Base structure
    T = [(0, 1)]
    lo, hi = 0, 1  # maintain current span [lo, hi]

    # Recursively refine k times
    for _ in range(k):
        # Normalize to start at 0 to simplify arithmetic and keep integers
        if lo != 0:
            T = [(l - lo, r - lo) for (l, r) in T]
            hi -= lo
            lo = 0

        delta = hi - lo  # current normalized span length (integer)

        S = []
        # 4 tiled copies of the current structure
        for start in tile_positions:
            offset = delta * start
            S.extend((l + offset, r + offset) for (l, r) in T)

        # 4 blocking intervals scaled to current span
        S.extend((delta * a, delta * b) for (a, b) in blocker_windows)

        # Update structure and its new span: it always lies within [2*delta, 16*delta]
        T = S
        lo, hi = 2 * delta, 16 * delta

    return T

# EVOLVE-BLOCK-END

def run_experiment(**kwargs):
  """Main called by evaluator"""
  return construct_intervals()