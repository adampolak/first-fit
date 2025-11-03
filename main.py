#!/usr/bin/env python3
from shinka.core import EvolutionRunner, EvolutionConfig
from shinka.database import DatabaseConfig
from shinka.launch import LocalJobConfig

job_config = LocalJobConfig(eval_program_path="evaluate_ratio.py")

strategy = "weighted"
if strategy == "uniform":
  # 1. Uniform from correct programs
  parent_config = dict(
    parent_selection_strategy="power_law",
    exploitation_alpha=0.0,
    exploitation_ratio=1.0,
  )
elif strategy == "hill_climbing":
  # 2. Hill Climbing (Always from the Best)
  parent_config = dict(
    parent_selection_strategy="power_law",
    exploitation_alpha=100.0,
    exploitation_ratio=1.0,
  )
elif strategy == "weighted":
  # 3. Weighted Prioritization
  parent_config = dict(
    parent_selection_strategy="weighted",
    parent_selection_lambda=10.0,
  )
elif strategy == "power_law":
  # 4. Power-Law Prioritization
  parent_config = dict(
    parent_selection_strategy="power_law",
    exploitation_alpha=1.0,
    exploitation_ratio=0.2,
  )
elif strategy == "power_law_high":
  # 4. Power-Law Prioritization
  parent_config = dict(
    parent_selection_strategy="power_law",
    exploitation_alpha=2.0,
    exploitation_ratio=0.2,
  )
elif strategy == "beam_search":
  # 5. Beam Search
  parent_config = dict(
    parent_selection_strategy="beam_search",
    num_beams=10,
  )

db_config = DatabaseConfig(
  db_path="evolution_db.sqlite",
  num_islands=2,
  archive_size=40,
  # Inspiration parameters
  elite_selection_ratio=0.3,
  num_archive_inspirations=4,
  num_top_k_inspirations=2,
  # Island migration parameters
  migration_interval=10,
  migration_rate=0.1,  # chance to migrate program to random island
  island_elitism=True,  # Island elite is protected from migration
  **parent_config,
)

search_task_sys_msg = """You are an expert mathematician and theoretical computer scientist specializing in online graph coloring algorithms, interval graphs, and computational geometry.

Try to construct a sequence of open intervals that maximizes the competitive ratio of the FirstFit algorithm. That is the number of colors used by FirstFit divided by the number of colors used by an optimal offline algorithm should be maximized.

Note that the offline optimum is equal to the clique number of the intersection graph of the intervals, which in turn is equal to the largest number of intervals that cover a single point.

The best known lower bound is 5. The best known upper bound is 8.

Try using recursive strategies.

Try getting inspired by this paper: https://arxiv.org/abs/1506.00192

For further inspiration try searching scientific literature on various variants of online interval coloring and related problems.

Heuristics that tend to raise FF without raising omega (suggestions)
- Keep a low omega "spine" (e.g., 10 <= omega <= 25) of long intervals, and inject waves of short intervals that pairwise avoid forming larger cliques but overlap many active colors so FF is forced upward.
- Use arrival‑order engineering: place blocker intervals to occupy small colors early; later intervals that overlap these blockers but mutually avoid forming a bigger clique push FF to assign new colors beyond 5*omega.
- Build towers: layered stacks of intervals with staggered starts/ends; then caps that overlap one piece from each tower to couple colors across layers.
- Periodically run the shrinker to keep n and span modest; smaller witnesses are easier to verify.

Try using a SAT solver or an LP solver in your code. The Python environment in which your code is run has installed:
- python-sat (with all optional dependencies),
- pulp (with all open source solver).

Be creative and try to find a new solution better than the best known lower bound."""

# Optimized with OpenAI Prompt Optimizer
search_task_sys_msg_optimized = """# Role and Objective
You are an expert in mathematics and theoretical computer science, specializing in online graph coloring algorithms, interval graphs, and computational geometry. Your objective is to construct a sequence of open intervals that maximizes the competitive ratio of the FirstFit algorithm—specifically, to maximize the ratio of colors used by FirstFit compared to the optimal offline solution.
Begin with a concise checklist (3–7 bullets) of your sub-tasks before starting substantive work; keep the items conceptual and high-level.
# Instructions
- Construct a sequence of open intervals that maximizes FirstFit's competitive ratio (number of colors used by FirstFit divided by the number used by the optimal offline algorithm).
- Recall: The offline optimum equals the clique number of the intersection graph, i.e., the largest number of intervals covering a single point (ω).
- Current best-known results:
- Lower bound: 5
- Upper bound: 8
- Apply recursive strategies where helpful.
- Draw inspiration from: https://arxiv.org/abs/1506.00192
- Additionally, consult the scientific literature for variants on online interval coloring and related problems.
## Heuristic Suggestions
- Maintain a low ω "spine" (e.g., 10≤ω≤25) of long intervals; inject waves of short intervals that pairwise avoid forming larger cliques but overlap many active colors, pushing FirstFit's color usage up.
- Engineer the arrival order: arrange for early "blocker" intervals to occupy low colors; subsequently, intervals overlapping these blockers but not each other can force FirstFit to introduce new colors above 5·ω.
- Build "towers": create layers of staggered intervals; then add capping intervals that each intersect one interval from every tower, coupling colors across layers.
- Use periodic "shrinking" strategies to keep both n (number of intervals) and total span modest, enabling easier verification and smaller witness constructions.
# Additional Guidelines
- Experiment with SAT or LP solvers in your implementation.
- Available packages in your Python environment:
- `python-sat` (with all optional dependencies)
- `pulp` (with all open-source solver backends)
- Strive for creativity and aim to construct a new example that improves on the current best-known lower bound.
# Context
- Focus on constructing sequences maximizing the competitive ratio of FirstFit for interval coloring.
- Reference literature such as https://arxiv.org/abs/1506.00192 for existing constructions and ideas.
- Incorporate heuristic and algorithmic suggestions into your approach, using available computational tools where suited.
# Post-action Validation
- After each substantive step (e.g., construction, algorithm proposal, code run, or solver invocation), briefly validate the outcome (in 1–2 lines) and decide on the next step or, if an issue arises, self-correct.
"""

evo_config = EvolutionConfig(
  task_sys_msg=search_task_sys_msg,
  # task_sys_msg=search_task_sys_msg_optimized,
  patch_types=["diff", "full", "cross"],
  patch_type_probs=[0.6, 0.3, 0.1],
  num_generations=400,
  max_parallel_jobs=5,
  max_patch_resamples=3,
  max_patch_attempts=3,
  job_type="local",
  language="python",
  llm_models=[
    # "gemini-2.5-pro",
    # "gemini-2.5-flash",
    # "bedrock/us.anthropic.claude-sonnet-4-20250514-v1:0",
    "o4-mini",
    "gpt-5",
    "gpt-5-mini",
    "gpt-5-nano",
    # "gpt-5-pro",
  ],
  llm_kwargs=dict(
    temperatures=[0.0, 0.5, 1.0],
    reasoning_efforts=["auto", "low", "medium", "high"],
    # reasoning_efforts=["high"],  # the only option possible for gpt-5-pro
    max_tokens=32768,
  ),
  meta_rec_interval=10,
  meta_llm_models=["gpt-5-nano"],
  meta_llm_kwargs=dict(temperatures=[0.0], max_tokens=16384),
  embedding_model="text-embedding-3-small",
  code_embed_sim_threshold=0.995,
  novelty_llm_models=["gpt-5-nano"],
  novelty_llm_kwargs=dict(temperatures=[0.0], max_tokens=16384),
  llm_dynamic_selection="ucb1",
  llm_dynamic_selection_kwargs=dict(exploration_coef=1.0),
  init_program_path="initial.py",
  results_dir="results",
)

if __name__ == "__main__":
  EvolutionRunner(
    evo_config=evo_config,
    job_config=job_config,
    db_config=db_config,
    verbose=True
  ).run()
