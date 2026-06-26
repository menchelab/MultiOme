"""Tests for the .mpx bundle round-trip and cross-validation."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from multiome_algo.crossval import retrieval_cv
from multiome_algo.groups import load_orphanet_groups
from multiome_core.bundle import read_bundle, write_bundle
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


def test_bundle_roundtrip(multiplex, tmp_path):
    bundle = write_bundle(multiplex, tmp_path / "slice.mpx")
    assert (bundle / "manifest.yaml").exists()

    loaded = read_bundle(bundle)
    assert loaded.layer_ids == multiplex.layer_ids
    for lid in multiplex.layer_ids:
        assert loaded[lid].n_edges == multiplex[lid].n_edges
        assert loaded[lid].n_nodes == multiplex[lid].n_nodes
        assert loaded[lid].scale == multiplex[lid].scale
    assert loaded.universe() == multiplex.universe()


def test_bundle_partial_load(multiplex, tmp_path):
    bundle = write_bundle(multiplex, tmp_path / "slice.mpx")
    loaded = read_bundle(bundle, layer_ids=["ppi"])
    assert loaded.layer_ids == ["ppi"]


def test_cv_informed_beats_or_matches_random(multiplex, groups):
    group = max(groups, key=len)
    rng = np.random.default_rng(0)
    res = retrieval_cv(
        multiplex, group, k_folds=3, weighting="informed",
        r=0.7, lcc_trials=50, rng=rng,
    )
    # informed propagation should retrieve held-out genes better than chance
    assert res.auroc > 0.5
    assert res.n_folds == 3
    assert 5 in res.top_k
