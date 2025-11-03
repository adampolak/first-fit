import argparse

from shinka.core import run_shinka_eval

MAX_NUM_INTERVALS = 10000

def simulate_first_fit(intervals: list[tuple[float, float]]) -> int:
  colors = []
  for (a_left, a_right) in intervals:
    used = set()
    for color, (b_left, b_right) in zip(colors, intervals):
      if not (a_right <= b_left or b_right <= a_left):
        used.add(color)
    color = 1
    while color in used:
      color += 1
    colors.append(color)
  return max(colors)

def simulate_offline_opt(intervals: list[tuple[float, float]]) -> int:
  events = []
  events += [(left, 1) for (left, _) in intervals]
  events += [(right, -1) for (_, right) in intervals]
  count = 0
  max_count = 0
  for x, v in sorted(events):
    count += v
    max_count = max(max_count, count)
  return max_count

def validate(intervals: list[tuple[float, float]]) -> tuple[bool, str]:
  if not isinstance(intervals, list):
    return False, "Not a list"
  if len(intervals) == 0:
    return False, "List of intervals cannot be empty"
  if len(intervals) > MAX_NUM_INTERVALS:
    return False, "Length %d exceeds the max limit of %d intervals" % (len(intervals), MAX_NUM_INTERVALS)
  for interval in intervals:
    if not isinstance(interval, tuple):
      return False, "Interval should be represented as a tuple"
    if len(interval) != 2:
      return False, "Not an interval, has more than 2 elements"
    left, right = interval
    if not (isinstance(left, float) or isinstance(left, int)):
      return False, "Left endpoint is neither float nor int"
    if not (isinstance(right, float) or isinstance(right, int)):
      return False, "Right endpoint is neither float nor int"
    if left >= right:
      return False, "Not an interval, the first element is not to the left of the right one"
  return True, "OK"

def compute_metrics(results: list[list[tuple[float, float]]]) -> dict:
  assert 1 == len(results)
  intervals = results[0]
  ALG = simulate_first_fit(intervals)
  OPT = simulate_offline_opt(intervals)
  competitive_ratio = ALG / OPT
  return {
    'combined_score': competitive_ratio,
    'public': {
      'num_intervals': len(intervals),
      'intervals': str(intervals)[:1000],
      'alg': ALG,
      'opt': OPT,
    }
  }

def main(program_path: str, results_dir: str):
  metrics, correct, error_msg = run_shinka_eval(
    program_path=program_path,
    results_dir=results_dir,
    experiment_fn_name="run_experiment",
    num_runs=1,
    validate_fn=validate,
    aggregate_metrics_fn=compute_metrics,
  )

if __name__ == "__main__":
  parser = argparse.ArgumentParser(
    description="Online interval coloring lower bound evaluator using shinka.eval"
  )
  parser.add_argument(
    "--program_path",
    type=str,
    default="initial.py",
    help="Path to program to evaluate (must contain 'run_experiment')",
  )
  parser.add_argument(
    "--results_dir",
    type=str,
    default="results",
    help="Dir to save results (metrics.json, correct.json)",
  )
  
  parsed_args = parser.parse_args()
  main(parsed_args.program_path, parsed_args.results_dir)
