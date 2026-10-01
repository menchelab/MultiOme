# MultiOme

Multiplex network medicine for **gene groups**. This is a Python implementation of the
method in Buphamalai et al., *Nature Communications* 2021
([10.1038/s41467-021-26674-1](https://www.nature.com/articles/s41467-021-26674-1)).

You supply network layers (protein interactions, tissue co-expression, pathways,
ontologies and so on) and a **gene group**: any gene set, such as a rare-disease gene set,
a clinical panel, a GO term or your own list. MultiOme then does three things:

1. Measures how **modular** the group is in each layer, using the largest connected
   component (LCC) z-score against random gene sets with BH correction.
2. Runs an **informed multiplex random walk with restart**. Only the layers that are
   significant for the group are walked, and the walker spends time in each layer in
   proportion to its z-score. The output is a ranking of candidate genes.
3. Evaluates retrieval with **10-fold cross-validation** (AUROC, top-k), paired across
   configurations: informed vs all layers vs PPI only.

It also makes the plots for each step: a modularity heatmap, layer relevance per group,
CV performance, and the candidate ranking with per-layer contributions.

By default the method behaves as in the paper, with modern options available as flags.
Everything that differs from the original R code, including flaws found in it, is listed
in **[docs/DEVIATIONS.md](docs/DEVIATIONS.md)**.

## Install

```bash
uv sync --extra dev        # or: pip install -e ".[dev]"
uv run pytest
```

This installs the `multiome` command. Python ≥ 3.11. networkx is optional and only
needed for `Layer.to_networkx()`: `pip install -e ".[networkx]"`. Re-running
`tests/r_parity/make_golden.R` requires R; the golden files are committed.

## Input formats

**Network: a folder with one edge-list file per layer** (`.tsv`, `.csv`, `.txt`,
optionally gzipped). The file stem becomes the layer id.

- The first two columns are gene A and gene B. Edges are undirected; self-loops, NAs and
  duplicates are dropped.
- The delimiter is detected automatically (tab, comma, semicolon or whitespace).
- A header row is detected automatically.
- To use edge weights, put them in the third column and pass `--edge-weights`
  (`weights=True` in Python). Without that, layers are treated as binary.
- Optionally add a `layers.tsv` (or `layers.csv`) with columns `layer_id`, `tags` and
  `description`. Tags are `;`-separated free text and are used to group layers in plots.
- Parsed layers are cached as sparse matrices in `~/.cache/multiome/layers`. Turn this
  off with `--no-cache` or `cache_dir=None`.

**Gene groups**, in one of three formats:

- **GMT**: `name<TAB>description<TAB>gene<TAB>gene…`
- **Long table**: two columns, `group` and `gene`.
- **Wide table**: one row per group, with a column of genes separated by `;`, `,` or
  `|` (for example the Orphanet table in `data/`). The id column, label column and
  genes column are detected, or you can give them with `id_col`, `label_col` and
  `genes_col`.

**Gene ids.** Layers and groups must use the same namespace. For HGNC symbols,
`normalize_symbols()` (CLI: `--normalize-ids`) maps previous symbols, aliases and case
variants to current HGNC symbols. It uses an HGNC table downloaded once to
`~/.cache/multiome/`, and reports anything it cannot resolve. `coverage_report()` /
`coverage.tsv` show how many genes of each group are found in each layer.

## Python quickstart

```python
from multiome_core import read_multiplex, read_gene_groups, coverage_report, filter_groups
from multiome_core.ids import normalize_symbols
from multiome_algo import modularity_table, significant_layers, layer_weights, informed_rwr
from multiome_algo import retrieval_cv
from multiome_algo.viz import (plot_modularity_heatmap, plot_layer_weights,
                               plot_candidate_ranking, plot_cv_performance)

mpx = read_multiplex("my_layers/")                 # folder of edge lists
groups = read_gene_groups("my_groups.gmt")
mpx, groups, id_report = normalize_symbols(mpx, groups)   # optional, HGNC symbols
print(coverage_report(mpx, groups))
groups = filter_groups(groups, min_size=20, max_size=2000)

# 1. layer relevance: one row per group x layer, with z, p, BH q, significant
table = modularity_table(mpx, groups, n_trials=1000)
plot_modularity_heatmap(table, out="modularity.png")

# 2. informed RWR for one group
g = groups[0]
weights = layer_weights(significant_layers(table, g.id), "z")   # paper's weighting
result = informed_rwr(mpx, g.genes, layer_weights=weights, r=0.7)
print(result.top(20))
plot_layer_weights(table, g.id, out="layers.png")
plot_candidate_ranking(result, top_n=30, out="ranking.png")

# 3. cross-validation: informed vs all layers vs ppi (if present), same folds
cv = retrieval_cv(mpx, groups, k_folds=10, protocol="train")
print(cv.summary())
plot_cv_performance(cv, out="cv.png")
```

The paper's data ships in `data/`:

```python
from multiome_core.legacy import read_paper_dataset
mpx, groups, id_report = read_paper_dataset()      # 46 layers, 28 Orphanet groups
```

Main options (see the docstrings):

- `modularity_table(null="degree", p_value="empirical")`
- `informed_rwr(coupling="present", combine="geomean")`
- `layer_weights(..., "uniform" | "softmax")`
- `retrieval_cv(configs=[CVConfig(...)], protocol="paper")`

## Command line

```bash
# layer relevance for every group x layer  -> modularity.tsv, modularity_heatmap.png
multiome modularity --network my_layers/ --groups my_groups.gmt --out results/

# rank candidates for one group, or for a seed list
multiome rank --network my_layers/ --groups my_groups.gmt --group MY_SET --out results/ \
              --modularity results/modularity.tsv
multiome rank --network my_layers/ --seeds seeds.txt --name my_seeds --out results/

# cross-validated retrieval -> cv_folds.tsv, cv_summary.tsv, cv_performance.png
multiome cv --network my_layers/ --groups my_groups.gmt --out results/
multiome cv ... --configs informed=significant all=all:uniform ppi=ppi

# the same options from YAML (keys = option names; `command:` picks the subcommand)
multiome run config.yaml
```

Every command also writes `coverage.tsv`, and `symbol_mapping.tsv` when symbols are
normalised. Options shared across commands:

- `--layers`, `--min-size` / `--max-size`, `--seed`
- `--trials`, `--null {uniform,degree}`, `--p-value {norm,empirical}`, `--alpha`
- `--weighting {z,uniform,softmax}`, `--coupling {paper,present}`,
  `--combine {mean,geomean,rank_geomean}`, `--restart`

Run `multiome <command> --help` for the full list.

**Note on `rank`.** Without `--modularity`, layer significance is computed for that one
group, so the BH family is smaller than in the paper. Pass a `modularity.tsv` from a run
over all groups to use the paper's selection.

## Reproducing the paper

```bash
multiome modularity --paper --out paper/
multiome cv --paper --protocol paper --out paper_cv/     # as published
multiome cv --paper --out paper_cv_train/                # default, no leakage
```

`--paper` loads the shipped 46 layers and Orphanet groups and normalises symbols to
current HGNC (turn this off with `--no-normalize-ids`). Exact published numbers cannot be
regenerated from the original code. Why, and how results compare, is in
[docs/DEVIATIONS.md](docs/DEVIATIONS.md).

## Layout

- `src/multiome_core`: data model (`Layer`, `Multiplex`, `GeneGroup`), file readers,
  HGNC normalisation, paper-data loaders.
- `src/multiome_algo`: LCC modularity, the informed multiplex RWR, cross-validation,
  plots and the CLI.
- `tests/r_parity`: the original R functions and golden outputs used for parity tests.
- `docs/modernization/`: historical design notes. Network generation and new data sources
  are parked as future work.

## Citation

Buphamalai P, Kokotovic T, Nagy V, Menche J. Network analysis reveals rare disease
signatures across multiple levels of biological organization. *Nat Commun* 12, 6306
(2021). https://doi.org/10.1038/s41467-021-26674-1

## License

See [LICENSE.md](LICENSE.md).
