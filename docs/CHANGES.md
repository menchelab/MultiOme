# What changes from the original implementation

This package re-implements Buphamalai et al., *Nat Commun* 2021 in Python. The original R
code is at commit `cce867d` (`functions/`, `source/`). The goal is that the published
analysis can be reproduced on any network and gene set, so the paper's method is the
default. Improvements are either new defaults where they don't change the method's
meaning, or opt-in options. This page lists:

- **(a)** clarifications of how the method works, correcting earlier summaries (old README,
  `docs/modernization/`, the first Python scaffold);
- **(b)** what is kept exactly as in the paper;
- **(c)** improved defaults;
- **(d)** options that address known limitations, while the paper's behaviour stays the
  default.

Parity with R is tested: `tests/test_r_parity.py` runs the original R functions
(`tests/r_parity/original_functions.R`) on a toy multiplex. The layer transition matrix
matches to rtol 1e-12, and per-layer walk probabilities match to rtol 1e-6.

---

## (a) Clarifications

| Earlier summary | How the method works |
|---|---|
| Layer weights = softmax over LCC z-scores | Raw z-scores of the *significant* layers feed `pmat_cal` (below). Softmax is available as an option (`weighting="softmax"`). |
| A layer is informative if z ≥ 1.645 | p = `pnorm(z, lower=F)`. BH correction over all group×layer pairs, q < 0.05, plus ≥ 10 group genes in the layer and an LCC ≥ 5 (`process_LCC_result.R`). |
| 45 layers | **46** layers in `data/network_edgelists/` (~21.4 M edges). |
| PPI = BioPlex + HuRI | `ppi.tsv` is **HIPPIE**. BioPlex/HuRI files in `raw_data.zip` serve the PPI-subset analyses. |
| A global inter-layer jump parameter (`delta`, MultiXrank-style) | Inter-layer moves come from the weight-derived matrix `pmat_cal`, followed by a global column normalisation. |
| No multiple-testing correction | BH (`p.adjust(method="BH")`) is applied. |
| 5-fold CV with pooled AUROC (first Python scaffold) | **10-fold** CV. AUROC is computed per fold and summarised as median/IQR per disease. |

## (b) Kept as in the paper

| Step | Paper (R, `cce867d`) | Python default | Options |
|---|---|---|---|
| LCC modularity | `LCC_functions.R`: observed LCC vs 1000 random sets drawn uniformly from the layer's nodes; z = (obs − mean)/sd (sd with ddof = 1) | `modularity_table(..., n_trials=1000, null="uniform")` | `null="degree"`, `p_value="empirical"` |
| Significance | `process_LCC_result.R`: `pnorm` p, BH over all pairs, q < 0.05, n ≥ 10, LCC ≥ 5 | `p_value="norm"`, `alpha=0.05`, `min_genes=10`, `min_lcc=5` | — |
| Layer selection | the group's significant layers | `significant_layers(table, group_id)` | `CVConfig(layers="all" \| [ids])` |
| Layer transition matrix | `RWR_transitional_matrix.R::pmat_cal`: `P[i,j] = min(1, w_i/w_j)/L`, diagonal = 1 − Σ off-diagonal, w = z | `pmat_paper`, `weighting="z"` | `"uniform"`, `"softmax"`, explicit `P=` |
| Supra matrix | block (i,i) = `P[i,i]·M_i`, block (i,j) = `P[i,j]·I`, then column-normalised; every gene has a copy in every layer | `SupraOperator(coupling="paper")`: matrix-free, numerically identical | `coupling="present"` |
| Walk | `RWR.R`: r = 0.7; p0 = seeds replicated over layer copies, normalised | `informed_rwr(r=0.7)` | seed weights `{gene: w}` |
| Score | `weighted_multiplex_propagation.R`: arithmetic mean over layers, seeds removed | `combine="mean"`, `remove_seeds=True` | `"geomean"`, `"rank_geomean"` |
| Cross-validation | `10fold_retrieval_comparison.R`: 10 folds, 20–2000-gene groups; informed vs all layers uniform vs PPI | `retrieval_cv(k_folds=10, min_size=20, max_size=2000)`, `paper_configs()` | custom `CVConfig`s, `top_k` |

## (c) Improved defaults

**Layer selection from training genes in CV (`protocol="train"`).** Layers and weights
are chosen from each fold's training genes, which minimises leakage from held-out genes
into the model. The paper selected them once from the full group; `--protocol paper`
reproduces that. On the shipped data the two give almost the same median AUROC
(0.906 vs 0.911; see below).

**Gene symbol normalisation (shipped data).** The shipped layers were built at different
times and use different HGNC versions:

- the GTEx co-expression layers use SEPTIN7 and MARCHF5;
- PPI, GO, HP, MP and co-essential use SEPT7 and MARCH5;
- the Orphanet table has case variants (C9ORF72).

`normalize_symbols()` maps every symbol to the current HGNC symbol, in this order:
current, previous, alias, then case-insensitive variants. Ambiguous hits are reported and
left unchanged. On the shipped data, 825 symbols are renamed and 160 stay unresolved, and
unmatched group memberships drop from 154 to 70. Every run writes `symbol_mapping.tsv` and
`coverage.tsv`. For your own data, normalisation is off unless you pass `--normalize-ids`.

**Edge-list cleaning.** Missing values, self-loops and duplicate undirected edges are
dropped; for duplicates the maximum weight is kept. In `ppi.tsv` this removes 4,345
self-loops, 1,074 duplicate pairs and 774 NA rows. Quoted `"NA"` tokens in
`reactome_copathway.tsv` are read as missing values.

**Other:**

| Topic | Now |
|---|---|
| z when sd = 0 | `NaN`: the layer is reported as not assessable |
| Convergence | L1 change < `tol=1e-10` (max 1000 iterations) |
| Reproducibility | every null and fold is seeded from `(seed, layer/group id)`, independent of processing order |
| Folds | built once from group ∩ multiplex universe and **shared** across configurations, so comparisons are paired |
| Candidates / positives | candidates = genes of the configuration's walked layers minus seeds. Held-out genes absent from those layers are counted in `n_held_out_missing`. |
| Data packaging | plain files (edge-list folder, optional `layers.tsv`, GMT/TSV gene sets) plus a disk cache; layers carry free-text `tags` |

## (d) Options for known limitations

The paper's choices remain the default, so results stay comparable to the publication.

**Degree-preserving null (`null="degree"`).** With the uniform null, random sets are drawn
uniformly from a layer's nodes. Disease genes tend to be well studied and have higher
degree, which raises their LCC under a uniform null. The degree-preserving null draws
random sets matched by degree. For the largest group:

| Layer | z, uniform null | z, degree null |
|---|---|---|
| ppi | 10.51 | 2.51 |
| coex_TST | 2.92 | −0.06 |

**Empirical p-values (`p_value="empirical"`).** The normal approximation `pnorm(z)` is
fast and works well for large LCCs. For small, discrete LCC distributions it can be
optimistic. In a toy test, a random 20-gene set had LCC 6 against a null with mean 2.86
and sd 0.93. That gives p_norm = 0.0003, while the empirical p is 0.02 (50 trials).
The empirical option uses (#rand ≥ obs + 1)/(trials + 1).

**Present-only coupling (`coupling="present"`).** In the paper's formulation every gene
has a state in every layer. A gene absent from a layer passes its probability straight
on to its other copies. With `coupling="present"`, only existing (gene, layer) pairs are
states, so the stationary layer distribution follows the weights more closely.

**Scope of the BH correction.** The paper corrected over all group×layer pairs at once,
and so does `modularity_table` over the groups you pass:

- Under `protocol="train"`, each fold's table covers all groups' training sets.
- `multiome rank` with a single group corrects over that group's layers only. Pass a
  full-run `--modularity modularity.tsv` to reproduce the paper's selection.

**HP/MP layers (open question).** The HPO and MPO gene-similarity layers may share
annotation sources with the Orphanet gene–disease table. If that matters for your
question, compare runs with `--layers` excluding HP/MP.

**Comparing with published numbers.** The published figures were produced over several
iterations of the analysis code, so expect close qualitative agreement rather than
identical values. Per-fold values are kept in `CVResult.folds`, and `summary()` reports
median/IQR across folds.

---

## Full-scale reproduction

These runs use all 46 shipped layers and the 26 Orphanet groups of 20–2000 genes, with
symbols normalised. Each is 10-fold CV with 1000 null trials and seed 0, run with
`multiome cv --paper [--protocol paper] [--null degree]` on a single laptop CPU core.

| Setting | Informed | All layers, uniform | PPI only | Informed > uniform (groups) | Mean layers used | Runtime |
|---|---|---|---|---|---|---|
| `protocol="train"`, uniform null (default) | 0.906 | 0.858 | 0.723 | 19 / 26 | 11.2 | 18 min |
| `protocol="paper"`, uniform null (as published) | 0.911 | 0.858 | 0.723 | 20 / 26 | 12.1 | 7 min |
| `protocol="train"`, degree null | 0.887 | 0.858 | 0.723 | 17 / 26 | 10.1 | 67 min |

Values in the first three columns are the median over groups of each group's median
fold AUROC.

- **The paper's main claim holds.** Informed walks beat the all-layer walk by about 0.05
  median AUROC, and beat PPI alone by about 0.18.
- **Leakage in the published protocol is real but small**, about +0.005 AUROC.
- **The degree-preserving null selects fewer layers and lowers the informed AUROC** to
  0.887. That is still above the uniform walk on AUROC. Top-100 recovery is the exception:
  0.150 against 0.172 for the uniform walk (0.225 under the default null).
- **Early retrieval (top-10) barely separates the configurations** (0.026–0.032).
  Gains appear mainly over the top 100 and in the overall ranking.

