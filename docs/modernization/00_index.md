# MultiOme Modernization — Research Reference Index

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
