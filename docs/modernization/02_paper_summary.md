# Paper Summary — Buphamalai et al., Nat Commun 2021

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
> removed from the code. See `../DEVIATIONS.md` and the top-level README for current
> behaviour.

> "Network analysis reveals rare disease signatures across multiple levels of biological organization"
> Pisanu Buphamalai, Tomislav Kokotovic, Vanja Nagy, Jörg Menche. Nat Commun (2021). DOI: 10.1038/s41467-021-26674-1.
> Reference doc for modernization. Source: full-text read, 2026-06-26.

## Core Thesis
A multiplex network framework to investigate rare genetic diseases across biological scales. A 46-layer network spanning genome→phenome shows that the "disease module" concept (from complex-disease research) applies to rare diseases, and that **differential modularity across layers** both quantifies pathobiological relevance and enables accurate gene-candidate prediction.

## Scientific Questions
1. How do single-gene defects impact different scales between genotype and clinical phenotype?
2. Which network data types are most relevant for a given rare disease mechanism, and how to integrate heterogeneous omics?
3. Does the disease-module concept extend from polygenic to monogenic rare diseases?
4. Does modularity extend beyond PPI to transcriptome/pathway/phenotype scales?
5. Can insights from well-studied rare diseases inform poorly-characterized ones?
6. Can network methods prioritize causal genes in undiagnosed patients?

## Data Sources (7 major DBs + patient cohorts)
- **Genome**: CRISPR co-essentiality, 276 cancer cell lines (Kim et al. 2019).
- **Transcriptome**: GTEx v7 → 38 tissue-specific co-expression nets + 1 pan-tissue core.
- **Proteome**: HIPPIE v2.2 (2019) PPI.
- **Pathway**: Reactome (2019) co-membership.
- **Function**: GO-BP and GO-MF (2018) semantic similarity.
- **Phenotype**: HPO (human) + MPO (mouse) (2018) phenotypic similarity.
- **Disease associations**: Orphanet — 3,953 genes / 3,771 rare-disease terms → aggregated into **26 rare genetic disease groups**.
- **Patient cohorts**: local intellectual-disability cohort; RD-Connect GPAP (131 solved cases); temporal-holdout (21 patients, genes discovered after network construction).

## Methodology
**Multiplex construction**: 46 layers, 6 scales, 20M+ relationships among 20,354 genes. Three inference methods:
- Bipartite mapping (annotation-based, e.g. shared pathway).
- Semantic similarity (GO, HPO, MPO).
- Correlation (co-expression).
Disparity filter extracts significant edges from dense weighted nets. Co-expression: edges in ≤5 tissues = tissue-specific; others = core.

**Disease module detection (LCC)**: map disease-group genes to each layer; significance via z-score vs random gene sets; BH correction; p < 0.05.

**Embedding/viz**: node2vec 2D embeddings; MDS for layer similarity.

**Informed multiplex propagation** (the headline method): RWR on supra-adjacency matrix; identity matrices link layers; **inter-layer transition probability ∝ disease-specific layer relevance (π_dm)** derived from LCC modularity z-score (informative if BH-adjusted `pnorm` p < 0.05, with ≥ 10 genes and LCC ≥ 5). Detailed-balance condition makes the walker favor informative layers. Restart r = 0.7. Final visiting probability averaged across layers → gene ranking.

**Assessment**: 10-fold CV, AUROC primary metric, DeLong's test for ROC comparison. Benchmarked vs single-layer and vs gene-level features (pathway counts, expression, literature counts, phenotypic similarity).

## Multiplex Structure
- Genome (1): co-essentiality.
- Transcriptome (39): 38 tissue co-expression + 1 core.
- Proteome (1): HIPPIE PPI.
- Pathway (1): Reactome.
- Function (2): GO-BP, GO-MF.
- Phenotype (2): HPO, MPO.
Genes replicated across layers; identity matrices link a gene to itself across layers; inter-layer transition weights are disease-specific (modularity-driven).

## Key Findings
1. **Disease-module concept generalizes** to grouped rare diseases (avg 339 genes/group); 93% of individual rare diseases have <5 genes — below the ~20-gene detection threshold — hence grouping.
2. **Differential modularity**: PPI, HPO/MPO, core transcription, GO-BP modular across many groups; tissue-specific co-expression reveals disease-specific signatures (cardiac → heart nets); neuro diseases modular across many tissues (syndromic).
3. **Modularity predicts informative layers** → informed propagation beats single-layer.
4. **Gene prediction**: median CV AUROC **0.89** (informed multiplex) vs 0.81 (PPI-only), 0.85 (best single layer), 0.87 (all-layers equal). Dropping HPO → 0.80 (biggest single hit); dropping all curated layers → 0.71.
5. **Patient prioritization** (131 ID patients, mean 401 candidate genes): AUROC **0.95**; causal gene in **top 5 for 48.9%** (64/131). Gene-level baselines: top-5 in only 3.1–8.4%; phenotypic-similarity baseline AUROC 0.87.
6. **Temporal holdout** (21 patients): AUROC 0.86 (vs 0.90) — not just confirmatory bias.
7. **Syndromicity/size hurt performance**: more genes → lower AUROC (ρ = −0.83); multi-tissue modularity → lower (ρ = −0.53).

## Stated Limitations
1. Broad/heterogeneous disease-group definitions reduce specificity (esp. syndromic, multi-organ).
2. Data completeness: ~20-gene lower bound; PPI performance driven by size not curation detail.
3. Potential confirmatory bias from curated PPI/phenotype (mitigated by temporal holdout).
4. Performance degrades with group size & syndromicity — tension between enough genes and specificity.

## Future Directions
1. Finer, mechanism-based disease groupings; high-resolution phenotyping; subtyping syndromes.
2. Extend to cancer / complex disease; use modularity as a data-curation criterion.
3. Add omics layers: single-cell transcriptomics (cell-type modules), proteomics, metabolomics, epigenomics, spatial.
4. Combine propagation with deep learning; automate disease grouping/module discovery.
5. Clinical translation: prospective validation, integration with variant interpretation, beyond ID.

## Validation (multi-pronged)
- **Network-level**: node randomization (1000 sets); MDS/edge-overlap for layer distinctness; co-expression essentiality check; tissue-specific signal extraction.
- **Gene retrieval**: 10-fold CV on 26 groups; baselines (PPI, best-single, all-equal, informative-only); ablations (curated vs HT PPI; drop HPO/GO/Reactome; drop all curated).
- **Patient cohort**: 131 confirmed cases as gold standard; realistic post-filter candidate lists; vs gene-level features; AUROC + top-k + DeLong.
- **Temporal holdout**: 21 post-2019-discovery genes — bias-minimizing.
- **Biological**: disease–tissue correspondence; phenotypic coherence in HPO/MPO; known disease relationships recapitulated.

## Modernization-relevant takeaways
Strong foundation: multiplex construction + LCC modularity + informed RWR + patient validation. Opportunities: update data (GTEx v8+, newer PPI, current ontologies), add single-cell/spatial layers, explore DL alternatives to RWR, broaden beyond ID, prospective clinical deployment.
