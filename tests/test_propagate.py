from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse as sp

from multiome_algo.propagate import (
    SupraOperator,
    combine_layers,
    informed_rwr,
    layer_weights,
    normalized_adjacency,
    pmat_paper,
)


def test_pmat_stationary_distribution_proportional_to_weights():
    w = np.array([2.5, 4.0, 1.8])
    P = pmat_paper(w)
    np.testing.assert_allclose(P.sum(axis=0), 1.0)
    np.testing.assert_allclose(P @ (w / w.sum()), w / w.sum())  # detailed balance
    with pytest.raises(ValueError):
        pmat_paper([1.0, -0.5])


def test_layer_weights():
    z = {"a": 2.0, "b": 4.0}
    assert layer_weights(z, "z") == z
    assert layer_weights(z, "uniform") == {"a": 1.0, "b": 1.0}
    sm = layer_weights(z, "softmax")
    assert sm["b"] == 1.0 and np.isclose(sm["a"], np.exp(-2))


def test_single_layer_equals_plain_rwr(planted):
    mpx, groups = planted
    single = mpx.subset(["A"])
    op = SupraOperator(single, ["A"])
    seeds = sorted(groups[0].genes)[:10]
    X = op.walk(seeds, r=0.7, tol=1e-14)
    M = normalized_adjacency(single, "A")
    p0 = op.restart_vector(seeds)[:, 0]
    exact = sp.linalg.spsolve(sp.identity(M.shape[0], format="csc") - 0.3 * M.tocsc(), 0.7 * p0)
    np.testing.assert_allclose(X[:, 0], exact, atol=1e-12)


@pytest.mark.parametrize("coupling", ["paper", "present"])
def test_walk_conserves_mass(planted, coupling):
    mpx, groups = planted
    op = SupraOperator(mpx, mpx.layer_ids, {"A": 3.0, "B": 2.0, "C": 1.0}, coupling=coupling)
    X = op.walk(groups[0].genes)
    assert np.isclose(X.sum(), 1.0)
    assert (X >= -1e-15).all()
    if coupling == "present":
        assert np.allclose(X[~op.present], 0.0)  # no mass on absent copies


def test_informed_rwr_output(planted):
    mpx, groups = planted
    seeds = sorted(groups[0].genes)[:20]
    held = set(groups[0].genes) - set(seeds)
    res = informed_rwr(mpx, seeds, layer_weights={"A": 5.0, "B": 4.0})
    t = res.table
    assert t.loc[t["is_seed"], "rank"].isna().all()
    assert res.seeds == frozenset(seeds)
    assert list(res.layer_probs.columns) == ["A", "B"]
    top = set(res.top(len(held))["gene"])
    assert len(top & held) >= len(held) // 2  # planted module is recovered
    # combine options
    X = res.layer_probs.to_numpy()
    for how in ("mean", "geomean", "rank_geomean"):
        assert combine_layers(X, how).shape == (X.shape[0],)
    with pytest.raises(ValueError):
        combine_layers(X, "bogus")


def test_restart_rejects_absent_seeds(planted):
    mpx, _ = planted
    op = SupraOperator(mpx, ["A"])
    with pytest.raises(ValueError):
        op.restart_vector(["NOT_A_GENE"])
    X0 = op.restart_vector({"G000": 3.0, "G001": 1.0})
    i, j = op.index.get_indexer(["G000", "G001"])
    assert np.isclose(X0[i].sum() / X0[j].sum(), 3.0)
