# What deviates from the original implementation

This package re-implements Buphamalai et al., *Nat Commun* 2021 in Python. The original R
code is at commit `cce867d` (`functions/`, `source/`). This page lists:

- **(a)** errors in earlier descriptions of the method (old README, `docs/modernization/`,
  the first Python scaffold);
- **(b)** behaviour kept faithful to the paper by default;
- **(c)** defaults that were changed, and why;
- **(d)** flaws in the original method that are kept as the default for reproducibility,
  with an option to avoid them.

Parity with R is tested: `tests/test_r_parity.py` runs the verbatim R functions
(`tests/r_parity/original_functions.R`) on a toy multiplex. The layer transition matrix
matches to rtol 1e-12, and per-layer RWR probabilities match to rtol 1e-6.

---

## (a) Corrections to earlier descriptions

| Earlier claim | What the original code actually does |
|---|---|
| Layer weights = softmax over LCC z-scores, `exp(z)/Σexp(z)` | **No softmax anywhere.** Raw z-scores of the *significant* layers feed `pmat_cal` (below). Softmax is now only an opt-in option (`weighting="softmax"`). |
| A layer is informative if z ≥ 1.645 | Significance is p = `pnorm(z, lower=F)` (normal approximation). BH correction is applied over all group×layer pairs, q < 0.05, plus ≥ 10 group genes in the layer and an LCC ≥ 5 (`process_LCC_result.R`). |
| 45 layers | **46** layers in `data/network_edgelists/` (~21.4 M edges). |
| PPI = BioPlex + HuRI | `ppi.tsv` is **HIPPIE**. The BioPlex/HuRI files in `raw_data.zip` serve only the PPI-subset analyses. |
| A global inter-layer jump parameter (`delta`, MultiXrank-style) | There is none. Inter-layer moves come from the weight-derived matrix `pmat_cal`, followed by a global column normalisation. |
| No multiple-testing correction | BH (`p.adjust(method="BH")`) is applied. |
| 5-fold CV with pooled AUROC (first Python scaffold) | **10-fold** CV. AUROC is computed per fold and summarised as median/IQR per disease. |

## (b) Faithful by default

| Step | Original (R, `cce867d`) | Python default | Options |
|---|---|---|---|
| LCC modularity | `LCC_functions.R`: observed LCC vs 1000 random sets drawn uniformly from the layer's nodes, z = (obs − mean)/sd (sd with ddof = 1) | `modularity_table(..., n_trials=1000, null="uniform")` | `null="degree"`, `p_value="empirical"` |
| Significance | `process_LCC_result.R`: `pnorm` p, BH over all pairs (applied before the LCC ≥ 5 filter), q < 0.05, n ≥ 10, LCC ≥ 5 | `p_value="norm"`, `alpha=0.05`, `min_genes=10`, `min_lcc=5` | — |
| Layer selection | only the significant layers of the group | `significant_layers(table, group_id)` | `CVConfig(layers="all" \| [ids])` |
| Layer transition matrix | `RWR_transitional_matrix.R::pmat_cal`: `P[i,j] = min(1, w_i/w_j)/L`, diagonal = 1 − Σ off-diagonal, w = raw z | `pmat_paper`, `weighting="z"` | `"uniform"`, `"softmax"`, an explicit `P=` |
| Supra matrix | block (i,i) = `P[i,i]·M_i`, block (i,j) = `P[i,j]·I`; the whole matrix is then column-normalised; every gene has a copy in every layer | `SupraOperator(coupling="paper")`, matrix-free and numerically identical | `coupling="present"` |
| Walk | `RWR.R`: r = 0.7; p0 = seeds replicated over layer copies, normalised | `informed_rwr(r=0.7)` | seed weights `{gene: w}` |
| Score | `weighted_multiplex_propagation.R`: arithmetic mean over layers (`avg`), seeds removed | `combine="mean"`, `remove_seeds=True` | `"geomean"`, `"rank_geomean"` |
| Cross-validation | `10fold_retrieval_comparison.R`: 10 folds, 20–2000-gene groups, informed vs all-layers-uniform vs PPI | `retrieval_cv(k_folds=10, min_size=20, max_size=2000)`, `paper_configs()` | custom `CVConfig`s, `top_k` |

## (c) Changed defaults

**CV protocol: `protocol="train"` (default) vs `"paper"`.**
In the original, layer selection and weights came from the *full* gene group, held-out
genes included. The held-out genes therefore influenced which layers the walk used, which
leaks label information and inflates AUROC. The default now recomputes the modularity
table from each fold's training genes. `--protocol paper` reproduces the published setup.

**Symbol normalisation (on by default for the shipped data only).**
The shipped layers mix symbol versions:

- GTEx co-expression layers use new names (SEPTIN7, MARCHF5).
- PPI, GO, HP, MP and co-essential use old ones (SEPT7, MARCH5).
- The Orphanet table has case variants (C9ORF72).

Un-normalised, one gene becomes two disconnected nodes. `normalize_symbols()` maps to
current HGNC symbols in this order: current, previous, alias, then case-insensitive
variants. Ambiguous hits are reported and left unchanged.

On the shipped data, 825 symbols are renamed and 160 stay unresolved. Group memberships
with no matching layer gene drop from 154 to 70. Every run writes `symbol_mapping.tsv`
and `coverage.tsv`. For your own data, normalisation is off unless you pass
`--normalize-ids`.

**Edge-list cleaning.**
Missing values, self-loops and duplicate undirected edges are dropped; for duplicates the
maximum weight is kept. In `ppi.tsv` this removes 4,345 self-loops, 1,074 duplicate pairs
and 774 NA rows. The R-written quoted `"NA"` tokens in `reactome_copathway.tsv` are
treated as missing rather than as a gene called `"NA"`. Nodes whose only edge was a
self-loop disappear from that layer.

**Other changes:**

| Topic | Original | Now |
|---|---|---|
| z when sd = 0 | R division gives Inf/NaN | `NaN`: the layer is not assessable and never significant |
| Convergence | ad-hoc heuristic in `RWR.R` | L1 change < `tol=1e-10` (max 1000 iterations); parity test agrees to 1e-6 |
| Randomness | `seed=144` is an argument but never passed to `set.seed`, so runs are not reproducible | every null and fold is seeded from `(seed, layer/group id)`, independent of processing order |
| Folds | built per configuration over different gene universes | built once from group ∩ multiplex universe, **paired** across configurations |
| Candidates / positives | per-script | candidates = genes of the configuration's walked layers minus seeds. Held-out genes absent from those layers are reported (`n_held_out_missing`), not scored. |
| Data packaging | — | plain files (edge-list folder, optional `layers.tsv`, GMT/TSV gene sets) plus a disk cache. The earlier bundle/manifest and fixed biological-scale enum were dropped; layers carry free-text `tags`. |
| Absent genes | pass-through copies (see d) | unchanged by default; `coupling="present"` option |

## (d) Known flaws kept as default (with options)

**Uniform null inflates z.** Random sets are drawn uniformly from a layer's nodes. Disease
genes are better studied and have higher degree, so they form large LCCs by chance. A
degree-preserving null (`null="degree"`, `--null degree`) changes z drastically. For the
largest group:

| Layer | z, uniform null | z, degree null |
|---|---|---|
| ppi | 10.51 | 2.51 |
| coex_TST | 2.92 | −0.06 |

**Normal approximation of p.** `pnorm(z)` assumes normally distributed LCC sizes. They are
discrete and skewed, so p-values far beyond what 1000 trials can support are produced.
Option: `p_value="empirical"`, which gives (#rand ≥ obs + 1)/(trials + 1).

The approximation is **anticonservative for small LCCs**. In a toy test, a *random*
20-gene set had LCC 6 against a null with mean 2.86 and sd 0.93 (discrete and
right-skewed). That gives z = 3.39 and p_norm = 0.0003, while the empirical p is 0.02
(50 trials). Random sets can therefore be called significant.

**Pass-through coupling.** A gene absent from layer *j* still has a state in layer *j*.
That state has no intra-layer edges, so after column normalisation all its mass jumps
straight back out. This bends the intended stationary layer distribution and enlarges
the state space. Option: `coupling="present"`, where only (gene, layer) pairs that exist
are states and jumps go only to layers containing the gene.

**The BH family depends on the call.** The paper corrected over all group×layer pairs at
once. The same rule here means:

- Under `protocol="train"`, the family is each fold's table, which covers all groups'
  training sets.
- `multiome rank` with a single group corrects over that group's layers only, so it is
  less conservative than the paper. Pass a full-run `--modularity modularity.tsv` to
  reproduce the paper's selection.

**Possible circularity in HP/MP (hypothesis, not established).** The HPO and MPO
gene-similarity layers may share provenance with the Orphanet gene–disease annotations
used as groups. That would inflate those layers' modularity and their retrieval
performance. It has not been verified. Compare runs with `--layers` excluding HP/MP if
it matters for your conclusions.

**Reported CV variability (Supplementary Fig. 7d).** Per-disease mean AUROCs were plotted
as if they were folds, which understates variability. Here `CVResult.folds` keeps real
per-fold values; `summary()` reports median/IQR across folds.

**The committed code cannot regenerate the cached results.** The shipped `.RDS` caches were
produced by mixed code versions with different gene universes per configuration. The
size-grid/loess null in `LCC_CV_retrieval.R` is not used by the published pipeline, and
that script does not run as committed. Expect qualitative, not exact, agreement with the
published numbers.

---

## Full-scale reproduction

<!-- RESULTS -->
