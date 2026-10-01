# MultiOme Repository — Technical Analysis

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

> Reference doc for the modernization effort. Snapshot of the existing 2021 codebase: what it does and how. Source: automated repo analysis, 2026-06-26.

## 1. Repository Structure

```
data/        Raw + processed datasets, network edgelists
functions/   Core R function libraries (algorithms)
source/      Main analysis scripts (preprocess, compute, cross-validate)
report/      RMarkdown reproducible reports
Explorer/    Interactive Shiny web app for results
cache/       Pre-computed heavy results (~2.5 GB, downloaded separately)
deploy/      Deployment config
```

Top-level: `README.md`, `Dockerfile` (for the Shiny app only), `Multiome.Rproj`.

## 2. Data Inputs

### Raw data (`data/raw_data.zip`)

| Dataset | File | Purpose | Vintage |
|---------|------|---------|---------|
| GTEx | `GTEx_..._v7_..._gene_median_tpm.gct` (14.5 MB) | Gene expression → co-expression networks | v7 (2016) |
| HPO | `HPO_phenotype_to_genes.txt` (72 MB) | HPO gene–phenotype associations | Jun 2021 |
| OrphaNet | processed → `orphaNet_disease_gene_association_with_roots.tsv` | Rare disease–gene associations | — |
| BioPlex | `BioPlex.tsv`, `Bioplex_293T-*`, `Bioplex_HCT116-*` | PPI | BioPlex 3.0 |
| HuRI | `HuRI.tsv` (1.7 MB) | Reference interactome PPI | — |
| Reactome | `ReactomePathways.gmt` (853 KB) | Pathway co-membership | Jan 2019 |
| PanelApp | `PanelApp_ID_v3.0_2019-12-10.tsv` | Clinical gene panels | Dec 2019 |
| OGEE | `OGEE_esential_genes_20190416.txt` (11 MB) | Essential genes | Apr 2019 |

**Note:** `phenotype_annotation.tab` is listed in the zip's `readme.txt` but is NOT in the archive and is NOT referenced anywhere in code. Its description ("cached HPO-gene association before 2018") is also wrong — it is the HPO→disease annotation file (now `phenotype.hpoa`). Dead/ghost entry.

### Processed network edgelists (`data/network_edgelists/`) — 46 layers

1. **Co-expression (38 tissue-specific + 1 core)** from GTEx, disparity-filtered (corr ≥ 0.75, p < 0.01), tissue-specificity filtered (edges in ≤5 tissues). Files `coex_*.tsv`, plus `coex_core.tsv`. ~66K–1M edges each.
2. **PPI**: `ppi.tsv` (385K edges) — HIPPIE (not BioPlex/HuRI; those raw files serve only the PPI-subset analyses).
3. **Functional**: `GOBP.tsv` (179K), `GOMF.tsv` (19K).
4. **Phenotypic**: `HP.tsv` (84K), `MP.tsv` (34K).
5. **Other**: `co-essential.tsv` (68K, CRISPR), `reactome_copathway.tsv`.

### Disease/patient data
- `table_disease_gene_assoc_orphanet_genetic.tsv` — 20–2000 genes per disease group.
- `patient_hpo_terms.csv`, `patient_gene_list.csv` — for validation.

## 3. Core Methodology & Algorithms

### A. LCC (Largest Connected Component) modularity — `functions/LCC_functions.R`
`LCC_randomisation_measure(graph, nodesets, minnode=10, trial=1000)`:
1. Induce subgraph of disease genes; measure LCC size.
2. Sample 1000 random gene sets of size N.
3. z = (LCC_obs − mean_rand) / sd_rand; BH-FDR correction.
4. Flags which layers are significantly modular per disease (p < 0.05).

### B. Multiplex Random Walk with Restart — `functions/RWR.R`, `functions/weighted_multiplex_propagation.R`
Core: `p_{t+1} = (1−r)·W·p_t + r·p_0`, r = 0.7, W = column-normalized adjacency.
Multiplex extension:
1. Build supra-adjacency matrix across layers.
2. Keep only layers significant for the disease; their raw LCC z-scores w give the layer transition matrix `pmat_cal`: `P[i,j] = min(1, w_i/w_j)/L`, diagonal = 1 − Σ off-diagonal (no softmax).
3. Propagate; aggregate across layers (arithmetic mean, geometric mean of probs, geometric mean of ranks).
4. Return ranked gene list.

### C. Network distance — `functions/distance_functions.R`
Closest (`d_c`) and shortest (`d_s`) distance measures, with 1000-sample randomization for z-score/p-value. Used for disease separation and module compactness.

### D. Cross-validation — `source/LCC_CV_test.R`, `source/LCC_CV_retrieval.R`
10-fold CV per disease: 9 folds as seeds → propagate → measure retrieval on held-out fold (AUC, precision@k). Compares multiplex (disease-specific layers) vs PPI-only vs all-layers-unweighted. Finding: disease-specific multiplex wins.

### E. Network overlap — `source/network_overlap_randomisation.R`
Jaccard edge-overlap between layer pairs, with label-shuffling permutation test. Quantifies complementarity vs redundancy.

## 4. Biological Questions / Scope
1. Do rare diseases show modularity across biological scales (PPI, co-expression, functional, phenotypic)?
2. Are network types complementary or redundant?
3. Do diseases show tissue-specific co-expression signatures?
4. Can multiplex propagation beat PPI-only for gene prioritization?
5. Can it prioritize causal genes from patient variant + phenotype data?

Scope: 3,771 OrphaNet rare-disease terms → ~26–50 groups (sparse annotations forced grouping). Focus: genetic/Mendelian rare diseases.

## 5. Tech Stack
- **R 3.6.3** (primary), Bash for batch orchestration.
- igraph, tidygraph/ggraph, tidyverse, Matrix (sparse), MASS, pROC.
- Viz: ggplot2, cowplot, patchwork, visNetwork, ggiraph.
- Reports: rmarkdown, knitr, DT. App: shiny.
- Docker for the Shiny app only. No `renv.lock`. `pbapply` for progress (no real parallelism).

## 6. Outputs
- Topology metrics, LCC results (disease×network matrix of z/p/FDR), CV rankings (multiplex vs baselines), patient prioritization ranks, overlap matrices — all `.RDS`.
- Figures via RMarkdown (network complementarity, modularity heatmaps, tissue contextualization, CV AUC, patient case studies).
- Shiny Explorer (3 panels: differential modularity, network landscape t-SNE, network-disease inspection). Deployed at menchelab.com/MultiOmeExplorer.
- Network format: 2-column TSV edgelists.

## 7. Pain Points / Dated Aspects
1. **No dependency management** (no renv.lock / requirements; only Dockerfile pins).
2. **Hardcoded absolute paths** (e.g. `~/Documents/projects/Multiome/` in `functions/main_localisation_measure_function_nodist.R`).
3. **R 3.6.3** (2020) vs current 4.4+.
4. **Dated data**: GTEx v7 (2016), BioPlex 3.0, ontologies 2018.
5. **Manual download** of 2.5 GB `cache/` from Google Drive; `raw_data.zip` committed to git.
6. **DB credentials in code** (`functions/ontology_similarity_functions.R`, MySQL readonly, host `menchelabdb.int.cemm.at`).
7. **No full-pipeline container** (only Shiny).
8. **Monolithic scripts** with copy-paste blocks (e.g. three near-identical blocks in `LCC_CV_test.R`).
9. **Inconsistent seeds**, silent `NA`-on-error, global vars.
10. **No tests**, sparse docstrings.

### Modernization opportunities flagged
renv / setup script / `targets` pipeline / CLI args / `future`+`furrr` parallelism / updated data sources / DVC for data / GitHub Actions / full containerization.

## Summary
Publication-quality multilayer network-medicine codebase. Core algorithms (LCC modularity + informed multiplex RWR) are sound. Weaknesses are typical academic-code issues: no env management, hardcoded paths, dated data, monolithic scripts, no tests.
