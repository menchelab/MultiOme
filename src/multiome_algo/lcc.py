"""Largest Connected Component (LCC) modularity for a gene group on a layer.

Ports `LCC_randomisation_measure` from the 2021 R code. Measures whether a gene
group's genes form an unexpectedly large connected subgraph within a network layer,
relative to random gene sets of the same size drawn from the layer's node universe.

The resulting z-score is the per-layer relevance signal that drives informed
multiplex propagation (see propagate.py).
"""

from __future__ import annotations

from dataclasses import dataclass

import networkx as nx
import numpy as np

from multiome_core.schema import GeneGroup, Layer


@dataclass
class LCCResult:
    layer_id: str
    group_id: str
    n_genes_in_layer: int  # group genes present in this layer
    lcc_size: int  # observed largest connected component size
    rand_mean: float
    rand_std: float
    z_score: float
    p_value: float  # empirical right-tailed p (fraction of randoms >= observed)


def _largest_cc_size(graph: nx.Graph, nodes: set[str]) -> int:
    """Size of the largest connected component in the subgraph induced by `nodes`."""
    sub = graph.subgraph(nodes)
    if sub.number_of_nodes() == 0:
        return 0
    return max((len(c) for c in nx.connected_components(sub)), default=0)


def lcc_modularity(
    layer: Layer,
    group: GeneGroup,
    n_trials: int = 1000,
    min_genes: int = 10,
    rng: np.random.Generator | None = None,
) -> LCCResult | None:
    """Compute LCC-size z-score for a gene group on a layer via degree-agnostic randomization.

    Random sets are sampled uniformly from the layer's node universe (matching the
    original implementation). Returns None if fewer than `min_genes` group genes are
    present in the layer (module too small to assess).
    """
    rng = rng or np.random.default_rng()
    universe = list(layer.graph.nodes())
    present = group.genes & set(universe)
    n = len(present)

    if n < min_genes:
        return None

    observed = _largest_cc_size(layer.graph, present)

    universe_arr = np.array(universe, dtype=object)
    rand_sizes = np.empty(n_trials, dtype=float)
    for i in range(n_trials):
        sample = set(rng.choice(universe_arr, size=n, replace=False))
        rand_sizes[i] = _largest_cc_size(layer.graph, sample)

    rand_mean = float(rand_sizes.mean())
    rand_std = float(rand_sizes.std(ddof=1))
    z = (observed - rand_mean) / rand_std if rand_std > 0 else 0.0
    p = float((rand_sizes >= observed).sum() + 1) / (n_trials + 1)  # add-one smoothing

    return LCCResult(
        layer_id=layer.id,
        group_id=group.id,
        n_genes_in_layer=n,
        lcc_size=observed,
        rand_mean=rand_mean,
        rand_std=rand_std,
        z_score=z,
        p_value=p,
    )
