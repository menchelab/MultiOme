"""End-to-end smoke test of the phase-1 vertical slice on real repo data.

legacy edge lists -> Multiplex -> GeneGroup -> LCC modularity -> informed RWR.
Uses a small subset of layers/groups so it runs fast.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from multiome_algo.groups import load_orphanet_groups
from multiome_algo.lcc import lcc_modularity
from multiome_algo.propagate import informed_rwr, softmax_weights
from multiome_core.legacy import read_multiplex

DATA = Path(__file__).resolve().parents[1] / "data"
EDGELISTS = DATA / "network_edgelists"
GROUPS_TSV = DATA / "table_disease_gene_assoc_orphanet_genetic.tsv"

SLICE_LAYERS = ["ppi", "coex_core", "HP"]


@pytest.fixture(scope="module")
def multiplex():
    return read_multiplex(EDGELISTS, layer_ids=SLICE_LAYERS, name="slice")


@pytest.fixture(scope="module")
def groups():
    return load_orphanet_groups(GROUPS_TSV)


def test_legacy_reader_loads_layers(multiplex):
    assert len(multiplex) == len(SLICE_LAYERS)
    for lid in SLICE_LAYERS:
        layer = multiplex[lid]
        assert layer.n_nodes > 0
        assert layer.n_edges > 0


def test_groups_load(groups):
    assert len(groups) > 0
    g = groups[0]
    assert len(g.genes) > 0
    assert g.label


def test_lcc_modularity_runs(multiplex, groups):
    group = max(groups, key=len)  # largest group -> genes present in layers
    res = lcc_modularity(
        multiplex["ppi"], group, n_trials=200, rng=np.random.default_rng(0)
    )
    assert res is not None
    assert res.lcc_size >= 1
    assert res.n_genes_in_layer >= 10


def test_informed_rwr_ranks_universe(multiplex, groups):
    group = max(groups, key=len)
    rng = np.random.default_rng(0)
    zs = {}
    for lid in multiplex.layer_ids:
        res = lcc_modularity(multiplex[lid], group, n_trials=100, rng=rng)
        zs[lid] = res.z_score if res else 0.0
    weights = softmax_weights(zs)
    assert abs(sum(weights.values()) - 1.0) < 1e-9

    scores = informed_rwr(multiplex, group.genes, layer_weights=weights, r=0.7)
    assert len(scores) == len(multiplex.universe())
    total = sum(scores.values())
    assert total > 0
    # seed genes present in universe should tend to score highly
    ranked = sorted(scores, key=scores.get, reverse=True)
    top = set(ranked[: max(50, len(group) // 2)])
    present_seeds = group.genes & multiplex.universe()
    assert len(top & present_seeds) > 0
