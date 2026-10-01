# Proposed Architecture — Phase 1 (Python, gene-group generalized)

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

> Draft for discussion, 2026-06-26. Goal: clean, tool-friendly Python reimplementation that stays scientifically close to Buphamalai et al. 2021, with current data layers, network-building scripts, and plots — generalized so the unit of analysis is a **gene group** (any gene set), not specifically a rare-disease group.

## TWO INDEPENDENT UNITS (the central boundary)
The system splits into two separately-shippable units joined only by a network interchange format:

- **Unit A — Network Generation** (`multiome-net`): data source → network/edgelist. Free to evolve: new GTEx versions, regulatory-network inference (beyond plain correlation), better ontology similarity, or entirely different network data sources. Output = standardized multiplex bundle.
- **Unit B — The Algorithm** (`multiome-algo`): consumes a multiplex bundle + gene groups → LCC modularity, informed multiplex RWR, cross-validation, plots. Knows nothing about how networks were built.

**Phase 1 priority = Unit B against the EXISTING edge lists** already in `data/network_edgelists/` (46 layers). This de-risks the algorithm reimplementation first (validatable against the paper's own networks), and lets Unit A be modernized in parallel/after. We also define a **better interchange format** and ship a converter from the legacy TSVs.

### Interchange contract: the multiplex bundle
A directory bundle (proposed extension `.mpx/` or a zip):
```
mybundle/
  manifest.yaml          # multiplex metadata
  layers/
    ppi.parquet          # edge table: source, target, weight (typed, compressed)
    coex_BRC.parquet
    ...
  nodes.parquet          # optional node table: id, namespace (e.g. HGNC), aliases
```
`manifest.yaml` per layer: `id, scale (genome/transcriptome/proteome/pathway/function/phenotype), directed, weighted, node_namespace, source, build_params, version, date`. This is the ONLY thing Unit B depends on. It is also trivially exportable to MultiXrank's per-layer-TSV + YAML config and to Cytoscape SIF.

Rationale for Parquet over 2-col TSV: typed columns, ~5-10x smaller, fast columnar reads for 1M-edge layers, carries weights + provenance without separate sidecars. Legacy TSV reader/writer kept for compatibility.

## Design principles
1. **Library-first, tool-friendly.** Every stage is an importable function with typed inputs/outputs and a thin CLI. No hardcoded paths, no notebooks-as-pipeline. Suitable to expose as tools later.
2. **Gene group as the core abstraction.** A `GeneGroup` = (id, label, set of genes, optional provenance). Rare-disease groups, PanelApp panels, GO terms, custom sets all produce `GeneGroup`s. The method never assumes "disease."
3. **Layers are pluggable.** A `Layer` = (id, scale, igraph/networkx graph, builder metadata). Adding a layer = registering a builder; nothing else changes.
4. **Data acquisition is reproducible & cached.** ToolUniverse / direct downloads behind a `sources/` adapter layer with on-disk cache + provenance manifest (version, date, query). No committed data blobs.
5. **Separation of concerns**: `sources` → `build` (networks) → `analyze` (LCC, RWR, CV) → `viz`. Each runnable independently.

## Proposed package layout (two units, shared schema)
```
multiome_core/                 # shared: the interchange contract (depended on by both units)
  schema.py                    # Layer, Multiplex (bundle), GeneGroup, results — typed
  bundle.py                    # read/write .mpx bundles (parquet) + manifest
  legacy.py                    # read existing data/network_edgelists/*.tsv -> Layer/bundle

multiome_net/                  # UNIT A — network generation (modernize later / in parallel)
  sources/                     # data acquisition adapters
    hpo.py orphanet.py gtex.py string.py depmap.py reactome.py go.py panelapp.py
    base.py                    # cache, provenance manifest, retry, ToolUniverse client
  build/                       # construction -> Layer
    coexpression.py            # GTEx -> (correlation | regulatory inference) -> disparity filter
    ppi.py                     # STRING -> confidence-thresholded
    similarity.py              # ontology semantic similarity (GO/HPO/MPO) -> gene-gene
    coessential.py             # DepMap -> co-essentiality
    pathway.py                 # Reactome co-membership
  assemble.py                  # layers -> .mpx bundle
  cli.py

multiome_algo/                 # UNIT B — the algorithm (PHASE 1 FOCUS)
  groups/                      # GeneGroup construction (orphanet, panelapp, go, custom)
  lcc.py                       # LCC modularity z-score (randomization)
  distance.py                  # closest/shortest distance
  propagate.py                 # informed multiplex RWR (supra-adjacency + per-group weights)
  crossval.py                  # k-fold retrieval CV, AUROC, top-k
  viz/                         # modularity heatmap, ROC, network landscape
  cli.py

config/                        # YAML run configs (declarative)
tests/                         # unit + tiny-fixture integration
```
Unit B depends only on `multiome_core` (the bundle), never on `multiome_net`. Unit A produces bundles. In phase 1, bundles come from `multiome_core.legacy` over the existing edge lists.

## Stage contracts (the "clean APIs")
- `sources.<x>.fetch(...) -> raw table + provenance` (cached).
- `build.<layer>.build(raw, params) -> Layer` (networkx graph + metadata).
- `build.registry.assemble(layers) -> Multiplex`.
- `groups.<x>.load(...) -> list[GeneGroup]`.
- `analyze.lcc.modularity(layer, group, trials) -> z, p, lcc_size`.
- `analyze.propagate.informed_rwr(multiplex, seeds, layer_weights, r) -> ranked Series`.
  - `layer_weights` default = raw z-scores of significant layers into `pmat_cal` (the paper's π_dm; earlier draft wrongly said softmax). Pluggable: uniform, single-layer, custom.
- `analyze.crossval.retrieval_cv(multiplex, group, k, weighting) -> per-fold ranks, AUROC, top-k`.
- `viz.*` consume the typed result objects.

## Key technical choices to confirm
- **Graph lib**: `networkx` (readable, MultiXrank-compatible) vs `igraph`-python (faster, closer to original R). Lean networkx for phase 1; swap hot paths to scipy.sparse if needed. **(open)**
- **RWR engine**: our own scipy.sparse supra-adjacency RWR (full control of per-group weighting) vs MultiXrank as solver. Recommendation: **build our own ~minimal RWR first** (it's small, keeps the LCC-weighting native), benchmark MultiXrank later. **(open)**
- **Semantic similarity** (GO/HPO/MPO): reimplement Resnik/Lin/Jiang-Conrath in Python, or pull precomputed. Decide per ontology. **(open)**
- **Config format**: YAML describing the multiplex (layers + params), sources, and group sets — so a "run" is fully declarative. **(proposed)**

## Phase-1 scope (deliverable) — Unit B against existing edge lists
A reproducible algorithm pipeline that, from config alone:
1. Loads the **existing** `data/network_edgelists/` (46 layers) via `multiome_core.legacy` into an `.mpx` bundle (also proves the interchange format + converter).
2. Loads gene groups (start with the repo's OrphaNet groups for comparability + at least one non-disease source — e.g. a PanelApp panel or GO term — to prove the gene-group generalization).
3. Computes LCC modularity per layer per group.
4. Runs informed multiplex RWR + k-fold retrieval CV, reports AUROC/top-k.
5. Produces core plots (modularity heatmap, ROC, layer landscape).
6. Smoke test on a tiny fixture + validates qualitatively against the paper's reported behavior on the same networks.

**Unit A (network generation) is phase 2+**: ToolUniverse/STRING/DepMap fetch, GTEx v11, regulatory-network inference as an alternative to correlation co-expression, modern ontology similarity — each emitting the same `.mpx` bundle so Unit B is unchanged.

## Validation plan
Because we deliberately allow different/better inputs and methods, target the paper's *qualitative* claims, not exact AUROCs: informed multiplex > PPI-only > single-layer; tissue-specific co-expression adds signal; HPO is the highest-impact layer. Phase 1, running on the paper's own edge lists, should land closest to the published numbers — a useful correctness check before Unit A changes the inputs. Track AUROC deltas as a sanity band.

## What is explicitly deferred
Heterogeneous KG node types, sub-group resolution, single-cell/spatial layers, DL propagation, the Shiny→web app port. Within Unit A: regulatory-network inference, alternative network sources (kept in the design as pluggable builders, not built in phase 1).
