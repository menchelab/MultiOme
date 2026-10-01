# MultiXrank — Tool Evaluation

> **Historical design note (pre-2026-10 rebuild).** Kept for context; the implemented
> package is narrower and differs in places. Factual errors found later:
> - The original R code uses **no softmax**: raw LCC z-scores of the *significant* layers
>   feed `pmat_cal` (`P[i,j] = min(1, z_i/z_j)/L`).
> - Significance is **not** z ≥ 1.645. It is a `pnorm` p-value, BH over all
>   group×layer pairs, q < 0.05, ≥ 10 genes and LCC ≥ 5.
> - There are **46** layers, not 45.
> - `ppi.tsv` is **HIPPIE**, not BioPlex + HuRI.
> - The original has no global `delta`-style inter-layer jump parameter, so MultiXrank's
>   parameters do not map onto it directly.
>
> **Parked as future work:** network generation (`multiome_net`), ToolUniverse data
> sourcing, and the `.mpx` bundle/manifest with a fixed scale enum. All of these were
> removed from the code. See `../CHANGES.md` and the top-level README for current
> behaviour.

> Python package for RWR on heterogeneous multilayer networks.
> Repo: https://github.com/anthbapt/multixrank · Docs: https://multixrank-doc.readthedocs.io/
> Paper: Baptista, González, Baudot. "Universal multilayer network exploration by random walk with restart." Commun Phys 5, 170 (2022). DOI: 10.1038/s42005-022-00937-9.
> Companion: https://github.com/anthbapt/multixrank-tools · Maintainer: Anthony Baptista (KCL).
> Reference doc, 2026-06-26.

## What It Is
Explores **heterogeneous multilayer (universal multilayer) networks** via Random Walk with Restart. Starts from seed nodes and traverses arbitrarily complex architectures (multiple network types, each with multiple layers and node types) to rank nodes by network proximity to the seeds.

## Algorithm — RWR-MH (Multiplex + Heterogeneous)
A walker from seeds, at each step:
- Moves to a neighbor in the same layer (intra-layer).
- Jumps to another layer of the same multiplex — param **delta**.
- Jumps to another multiplex via bipartite edges — param **lambda**.
- Restarts from seeds — param **r**.
Iterates to steady state (residue < 1e-10). Global transition matrix integrates per-multiplex supra-adjacency matrices, bipartite matrices, and normalization. Final stationary distribution = node scores.

## Input Data Model
Supports **multiplex** (same nodes, different edge sets across layers) and **heterogeneous** (different multiplexes of different node types linked by bipartite edges). Directed/undirected, weighted/unweighted (graph_type codes 00/01/10/11).

Config = YAML + TSV edge lists:
```yaml
multiplex:
    1:
        layers:
            - multiplex/1/FR26.tsv
            - multiplex/1/FR3.tsv
    2:
        layers:
            - multiplex/2/UK15.tsv
bipartite:
    bipartite/1_2.tsv:
        source: 1
        target: 2
seed:
    seeds.txt
```
Edge lists: 2-col TSV (unweighted) or 3-col (weighted). Seeds: one ID per line.

## Outputs
Per-multiplex ranking TSV: `multiplex | node | score` (higher = closer to seeds). Optional SIF export of top-N subnetwork for Cytoscape.

## API / Usage
```python
import multixrank
multixrank.Example().write(path="airport")          # optional demo data
obj = multixrank.Multixrank(config="airport/config_minimal.yml", wdir="airport")
df = obj.random_walk_rank()
obj.write_ranking(df, path="output_airport")
obj.to_sif(df, path="output_airport/top3.sif", top=3)
```
Key params: **r** (restart), **delta** (intra-multiplex layer jump), **lambda** (inter-multiplex jump matrix), **eta** (restart dist across multiplexes), **tau** (restart dist across layers within a multiplex).

## Dependencies & Maturity
networkx, scipy, numpy, pandas, pyyaml. PyPI v0.3. Repo created Jul 2021; last commit Dec 2024; 20★/10 forks; MIT; Python 3.12; GitHub Actions CI; ReadTheDocs; peer-reviewed (Commun Phys 2022). Actively maintained, small community. Python/NumPy/SciPy — fine for moderate scale, not optimized for millions of nodes.

## Use Cases
Disease-gene prioritization, multi-omic integration, link prediction, leave-one-out CV. Biomedical KGs connecting genes/proteins/diseases/drugs/phenotypes. Companion `multixrank-tools` adds LOO CV, link-prediction benchmarks, parallelized parameter-space exploration.

## Relevance to MultiOme
**Strong fit** for the core method: MultiOme's hand-rolled `weighted_multiplex_propagation.R` is essentially an informed multiplex RWR — exactly MultiXrank's domain. Replacing the bespoke R supra-adjacency code with MultiXrank would:
- Remove custom matrix-assembly/normalization code.
- Add native heterogeneous support (could model genes AND diseases AND phenotypes as separate multiplexes joined by bipartite edges, rather than baking phenotype similarity into gene-gene layers).
- Provide tunable, documented, cited parameters.

**Caveats / gaps vs the paper's method:**
- MultiOme's key innovation is **disease-specific layer weighting from LCC modularity z-scores** (the π_dm weighting via `pmat_cal`, not softmax). MultiXrank's `delta`/`tau`/`eta`/`lambda` are global params, not per-disease learned weights. To reproduce the paper's "informed" propagation, layer weights would need to be injected per disease (set tau/eta from LCC z-scores per run) — feasible but requires wrapping MultiXrank in a per-disease loop that recomputes params.
- Scale: 46 layers × 20K nodes per disease group × many diseases × CV folds — verify performance is acceptable in Python vs the cluster-precomputed R.

**Verdict**: Best candidate to replace the custom propagation engine, *if* we can map LCC-derived layer relevance onto its parameters. Worth a prototype. Alternatively, keep LCC weighting logic ours and use MultiXrank only as the RWR solver.

## When NOT Needed
Single network → `networkx.pagerank` / scikit-network. Homogeneous multiplex with no cross-type nodes → simpler RWR. Need embeddings for downstream ML → node2vec/KG embeddings. Massive graphs → C++/GPU implementations.
