# ToolUniverse — Data-Access Platform Evaluation

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

> Unified API over 1000+ scientific databases. MCP server connected in this session.
> Reference doc, 2026-06-26. Tool counts/versions as reported by the live server.

## What It Is
A scientific tooling/data aggregation platform exposing **2,524 tools across 555 categories** through one API. Integrates 100+ biomedical/genomic/chemical/clinical/literature databases. Most need no API key. Standardized query patterns, cross-referencing built in.

Integrations include: OpenTargets, ChEMBL, Orphanet, HPO, GTEx, STRING, Reactome, KEGG, ClinVar, gnomAD, PubMed/EuropePMC, PDB, UniProt, Ensembl, Monarch, DepMap, PanelApp.

## Relevant Data-Access Tools (vs MultiOme's pipeline)

### HPO — 6 tools
`HPO_search_terms`, `HPO_get_term`, `HPO_get_genes_by_phenotype` (phenotype→genes), `HPO_get_diseases_by_phenotype`, `HPO_get_disease_annotations`. Returns gene symbols, NCBI IDs, disease IDs (ORPHA/OMIM/DECIPHER), MONDO xrefs.

### Orphanet — 10 tools
`Orphanet_search_diseases`, `Orphanet_get_disease`, `Orphanet_get_genes` (disease→gene, with causative/modifier/susceptibility type), `Orphanet_get_phenotypes` (disease→HPO + frequency), `Orphanet_get_gene_diseases`, `Orphanet_get_epidemiology`, `Orphanet_get_natural_history` (onset, inheritance), `Orphanet_get_icd_mapping`.

### GTEx — 16 tools (gtex + gtex_v2)
`GTEx_get_median_gene_expression` (median TPM, 54 tissues, **Adult GTEx V11, Jan 2026** — far newer than the repo's v7), `GTEx_get_expression_summary`, `GTEx_get_top_expressed_genes`, `GTEx_get_median_transcript_expression`, `GTEx_query_eqtl`.

### PPI
- **STRING — 6 tools**: `STRING_get_network` (edges + confidence by evidence type), `STRING_get_interaction_partners`, `STRING_functional_enrichment`, `STRING_get_functional_annotations`, `STRING_map_identifiers`, `STRING_ppi_enrichment`.
- **IntAct/Reactome**: `EBIProteins_get_interactions`, `ReactomeInteractors_get_protein_interactors`.
- **Others**: `NDEx_*` (published nets), `PathwayCommons_get_neighborhood` (22 DBs), `OmniPath_get_signaling_interactions`, `iPTMnet_get_ptm_ppi`.
- **BioPlex and HuRI: no dedicated tools.**

### Pathways
- **Reactome — 20 tools**: `Reactome_map_uniprot_to_pathways`, `Reactome_get_pathway`, `Reactome_get_pathway_hierarchy`, `Reactome_list_top_pathways`, `ReactomeAnalysis_species_comparison`.
- **KEGG — 12 tools**, **WikiPathways — 5**, **PathwayCommons — 4**.

### Gene Ontology — 5 tools
`GO_search_terms`, `GOAPI_get_genes_by_function`, `QuickGO_annotations_by_gene`.

### Essentiality
`OpenTargets_get_target_depmap_essentiality` (DepMap Chronos), `DepMap_get_gene_dependencies`. **OGEE: no dedicated tool** (DepMap is the closest analog).

### PanelApp — 3 tools
`PanelApp_search_panels`, `PanelApp_get_panel` (genes + confidence 3/2/1 = green/amber/red), `PanelApp_search_genes`.

### Disease–gene aggregators (new capability vs 2021)
`gather_disease_profile` (Orphanet+OMIM+DisGeNET+OpenTargets+OLS), `gather_gene_disease_associations` (cross-ref + concordance), OpenTargets (61 tools), Monarch (13), GenCC (3), ClinGen (8), Gene2Phenotype.

### Variant/clinical (new)
ClinVar (4), gnomAD (11), `annotate_variant_multi_source`, `MyVariant_query_variants`, Ensembl VEP / GenomeNexus / OpenCRAVAT.

## Coverage Summary vs MultiOme Inputs
| Source | ToolUniverse | Notes |
|--------|--------------|-------|
| HPO | ✓ (6) | full |
| OrphaNet | ✓ (10) | full |
| GTEx | ✓ (16) | **v11 (2026)** vs repo v7 |
| Reactome | ✓ (20) | full |
| PanelApp | ✓ (3) | full |
| GO | ✓ (5) | full |
| PPI (BioPlex) | ⚠ | no dedicated tool; STRING/IntAct overlap; or NDEx if uploaded |
| PPI (HuRI) | ⚠ | no dedicated tool; same workarounds |
| OGEE essential | ⚠ | use DepMap CRISPR essentiality instead |
| Co-essentiality (CRISPR) | ✓ | via DepMap |

## API Pattern
1. **Discover**: `find_tools(query="disease gene associations")` (semantic) or `list_tools(mode="categories")` / `grep_tools(pattern=..., field="description")`.
2. **Inspect**: `get_tool_info(tool_names=[...], detail_level="full")`.
3. **Execute**: `execute_tool(tool_name=..., arguments={...})`.

Example chain: `HPO_search_terms("seizure")` → `HPO_get_genes_by_phenotype("HP:0001250")` → `GTEx_get_median_gene_expression("SCN1A")` → `STRING_get_network("SCN1A")` → `Reactome_map_uniprot_to_pathways(...)`.

## Relevance to MultiOme
**Replaces most manual downloads + the `raw_data.zip` / Google-Drive-cache pattern** for HPO, OrphaNet, GTEx, Reactome, PanelApp, GO — with newer data. Enables a reproducible, scriptable data-acquisition layer (no committed binaries).

**Gaps to source elsewhere**: BioPlex, HuRI (use STRING/IntAct, or download BioPlex/HuRI directly), OGEE (use DepMap). Decide whether to switch the PPI layer to STRING (changes results, needs revalidation) or keep BioPlex+HuRI from primary sources.

**New opportunities**: disease-gene aggregators (GenCC/ClinGen/Monarch/OpenTargets) and variant tools (ClinVar/gnomAD/VEP) could modernize both the disease-gene ground truth and the patient variant-filtering step that the paper did manually.
