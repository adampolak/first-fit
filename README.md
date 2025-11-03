# Lower bound for FirstFit for interval coloring

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
```
