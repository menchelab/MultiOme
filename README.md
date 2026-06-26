# MultiOme

Multiplex network medicine for **gene groups** — a Python reimplementation and
generalization of the framework in Buphamalai et al., *Nature Communications* 2021
([10.1038/s41467-021-26674-1](https://www.nature.com/articles/s41467-021-26674-1)).

Given a set of network layers spanning multiple levels of biological organization
(protein interactions, tissue co-expression, pathways, ontologies, …) and a **gene
group** (any gene set — a rare-disease gene set, a clinical panel, a GO term, or a
custom list), it:

1. measures how **modular** the group is within each layer (LCC z-score), then
2. runs **informed multiplex random walk with restart** — propagation weighted by each
   layer's per-group relevance — to rank candidate genes.

## Architecture: two independent units

The system is split along one boundary, joined only by a network interchange format:

- **`multiome_net`** — *network generation* (phase 2+). Fetches data (ToolUniverse,
  STRING, DepMap, GTEx, …) and builds layers, emitting a standardized multiplex bundle.
  Free to evolve: new data versions, regulatory-network inference, new sources.
- **`multiome_algo`** — *the algorithm* (phase 1, current focus). Consumes a multiplex
  bundle + gene groups → LCC modularity, informed RWR, cross-validation, plots. Knows
  nothing about how the networks were built.
- **`multiome_core`** — the shared interchange contract: `Layer`, `Multiplex`,
  `GeneGroup`, bundle I/O, and a reader for the original 2021 TSV edge lists.

See [`docs/modernization/`](docs/modernization/) for the full design and the background
research (paper summary, repo analysis, MultiXrank and ToolUniverse evaluations).

## Status

Phase 1 runs the algorithm against the **original 46-layer edge lists** in
`data/network_edgelists/` (the paper's own networks) — the configuration where results
should land closest to the published numbers, used as a correctness check before
`multiome_net` modernizes the inputs.

## Setup

```bash
uv sync --extra dev
uv run pytest
```

## Quick start

```python
from multiome_core.legacy import read_multiplex
from multiome_algo.groups import load_orphanet_groups
from multiome_algo.lcc import lcc_modularity
from multiome_algo.propagate import informed_rwr, softmax_weights

mpx = read_multiplex("data/network_edgelists", layer_ids=["ppi", "coex_core", "HP"])
groups = load_orphanet_groups("data/table_disease_gene_assoc_orphanet_genetic.tsv")
group = groups[0]

zs = {lid: (r.z_score if (r := lcc_modularity(mpx[lid], group, n_trials=200)) else 0.0)
      for lid in mpx.layer_ids}
scores = informed_rwr(mpx, group.genes, layer_weights=softmax_weights(zs))
```

## License

See [LICENSE.md](LICENSE.md).
