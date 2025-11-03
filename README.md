# Lower bound for FirstFit for interval coloring

This is an experiment trying to use ShinkaEvolve in order to obtain a better lower bound for the competitive ratio of FirstFit for online interval coloring. The best known lower is 5 (https://arxiv.org/abs/1506.00192).

The `initial.py` solution in this repo produces an instance that proves 7/3=2.333..., and the best that ShinkaEvolve found (after burning $250 in OpenAI API calls) is 22/8=2.75.

## Quick start

```bash
# Install uv if you haven't already
curl -LsSf https://astral.sh/uv/install.sh | sh

# Create environment and install Shinka's dependencies
uv venv --python 3.11
source .venv/bin/activate
uv pip install -e 3rdparty/ShinkaEvolve

# Install LP and SAT solvers
uv pip install python-sat[aiger,approxmc,cryptosat,pblib]
uv pip install pulp[open_py]

# Set api key and url
export OPENAI_API_KEY=sk...
export OPENAI_BASE_URL=https://eu.api.openai.com/v1

# Run Shinka
python3 main.py
```

## Runs history

The `runs` folder contains logs from five most successful runs. Some files changed their names along the way, the code was cleaned up, etc., so they can be somewhat inconsistent with what the current code would produce.

- `runs/results` -- baseline, optimizes competitive ratio, run for ~200 iterations, the best solution is 13/5=2.6, found after 1 iteration;
– `runs/results_pro` -- same as above but with the gpt-5-pro model added to the list of models (uncomment  lines 89, 93, 94 in `main.py` to get it), run for ~50 iterations, results the same as previously;
– `runs/results_pro_prompt` -- same as above but with prompt optimized with OpenAI Prompt Optimizer (uncomment line 67 in `main.py` to get it), run for ~40 iterations, results still the same;
- `runs/results_numcolors` -- variant where instead of directly maximizing the competitive ratio we require opt<=10 and we maximize the number of colors used by FirstFit (uncomment lines 7 and 68 in `main.py`), run for 400 iterations, the best solution is 22/8=2.75, found after 91 iterations;


