# MultiOme Modernization — Research Reference Index

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

Background research compiled 2026-06-26 to inform modernizing the 2021 MultiOme codebase. These are reference docs, not the plan itself.

## Documents
1. [01_repo_analysis.md](01_repo_analysis.md) — What the existing codebase does and how (structure, data, algorithms, tech stack, pain points).
2. [02_paper_summary.md](02_paper_summary.md) — Buphamalai et al., Nat Commun 2021: questions, methods, multiplex structure, findings, limitations, stated future directions.
3. [03_multixrank.md](03_multixrank.md) — MultiXrank (RWR-MH) tool evaluation; fit/caveats as a replacement for the custom propagation engine.
4. [04_tooluniverse.md](04_tooluniverse.md) — ToolUniverse data-access platform; coverage of the repo's data sources and gaps.

## One-line synthesis
The paper's method = informed multiplex RWR with **disease-specific layer weights from LCC modularity z-scores**. MultiXrank can serve as the RWR engine but needs per-disease parameterization to reproduce the "informed" weighting. ToolUniverse can replace most manual data downloads (HPO, OrphaNet, GTEx v11, Reactome, PanelApp, GO) with current data; BioPlex/HuRI/OGEE need alternative sourcing (STRING/IntAct/DepMap or direct download).

## Decisions (2026-06-26)
- **Python end-to-end.** Move off R.
- **Phase 1 = clean foundation**: well-defined Python APIs usable as tools, current data layers, network-building scripts, plotting. Stay scientifically close to the original informed-multiplex-RWR method.
- **Generalize**: unit of analysis is a generic **gene group / gene set**, not specifically a rare-disease group. Rare disease becomes one use case.
- **Data via ToolUniverse** (HPO, OrphaNet, GTEx v11, Reactome, PanelApp, GO). **Switch PPI → STRING, essentiality → DepMap** (no ToolUniverse tools for BioPlex/HuRI/OGEE). Accept divergence from paper; revalidate.
- **Keep our LCC layer-weighting logic** (the novelty); MultiXrank only as a candidate RWR solver.

## Parked for later (keep in mind, not phase 1)
- Heterogeneous KG reframing (genes + diseases + phenotypes as distinct node types).
- Finer-than-group resolution.
- Single-cell / spatial / proteomics layers.
- DL alternatives to RWR.

See [05_architecture.md](05_architecture.md) for the proposed phase-1 design.
