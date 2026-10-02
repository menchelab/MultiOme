from __future__ import annotations

import networkx as nx
import numpy as np
import pandas as pd
import pytest

from multiome_algo.lcc import (
    LCCNull,
    bh_adjust,
    degree_bins,
    lcc_modularity,
    lcc_size,
    modularity_table,
    significant_layers,
)


def test_lcc_size_matches_networkx(planted):
    mpx, _ = planted
    layer = mpx["C"]
    G = layer.to_networkx()
    rng = np.random.default_rng(1)
    for _ in range(20):
        pos = rng.choice(layer.n_nodes, size=40, replace=False)
        sub = G.subgraph(layer.genes[pos])
        expected = max(len(c) for c in nx.connected_components(sub))
        assert lcc_size(layer.adj, pos) == expected


def test_bh_matches_r_p_adjust():
    # p.adjust(c(0.001, 0.01, 0.04, 0.2, 0.5, NA), "BH")
    p = np.array([0.001, 0.01, 0.04, 0.2, 0.5, np.nan])
    q = bh_adjust(p)
    np.testing.assert_allclose(q[:5], [0.005, 0.025, 0.04 * 5 / 3, 0.25, 0.5])
    assert np.isnan(q[5])
    # monotonicity step: p.adjust(c(0.01, 0.02, 0.03, 0.04, 0.05), "BH") == 0.05
    np.testing.assert_allclose(bh_adjust(np.array([0.01, 0.02, 0.03, 0.04, 0.05])), 0.05)


def test_planted_module_detected(planted):
    mpx, groups = planted
    tab = modularity_table(mpx, groups, n_trials=200)
    sig = significant_layers(tab, "module")
    assert set(sig) == {"A", "B"}
    assert all(z > 3 for z in sig.values())
    assert not tab.loc[tab["group_id"].str.startswith("rand"), "significant"].any()
    assert {"p_value", "q_value", "significant", "z_score"} <= set(tab.columns)


def test_null_is_order_independent_and_seeded(planted):
    mpx, groups = planted
    t1 = modularity_table(mpx, groups, n_trials=100, seed=3)
    t2 = modularity_table(mpx, groups[::-1], n_trials=100, seed=3)
    key = ["group_id", "layer_id"]
    a, b = t1.sort_values(key).reset_index(drop=True), t2.sort_values(key).reset_index(drop=True)
    np.testing.assert_allclose(a["z_score"], b["z_score"])
    t3 = modularity_table(mpx, groups, n_trials=100, seed=4)
    assert not np.allclose(t1["rand_mean"], t3["rand_mean"])


def test_min_genes_and_degree_null(planted):
    mpx, groups = planted
    assert lcc_modularity(mpx["A"], ["G000", "G001"], n_trials=10) is None
    res = lcc_modularity(mpx["A"], groups[0], n_trials=100, null="degree")
    assert res.z_score > 2
    bins = degree_bins(mpx["A"].degree, min_bin_size=50)
    counts = np.bincount(bins)
    assert counts.min() >= 50


def test_null_cache_reuse(planted):
    mpx, groups = planted
    null = LCCNull(mpx["A"], n_trials=50)
    pos = mpx["A"].positions(groups[0].genes)
    assert null.samples(pos) is null.samples(pos[::-1])


def test_invalid_args(planted):
    mpx, groups = planted
    with pytest.raises(ValueError):
        modularity_table(mpx, groups, p_value="bogus")
    with pytest.raises(ValueError):
        LCCNull(mpx["A"], null="bogus")


def test_heatmap_column_order(planted):
    from multiome_algo.viz.plots import _cluster_order, plot_modularity_heatmap

    mpx, groups = planted
    tab = modularity_table(mpx, groups, n_trials=20, seed=0)
    tags = {"A": "x", "B": "x", "C": "y"}
    fig = plot_modularity_heatmap(tab, layer_tags=tags)
    ax = fig.axes[0]
    assert [t.get_text() for t in ax.get_xticklabels()][-1] == "C"  # grouped by tag
    assert {"x", "y"} <= {t.get_text() for t in ax.texts}  # group labels drawn
    fig = plot_modularity_heatmap(tab, order=["C", "B", "A"])
    assert [t.get_text() for t in fig.axes[0].get_xticklabels()] == ["C", "B", "A"]
    with pytest.raises(ValueError):
        plot_modularity_heatmap(tab, order="bogus")
    # profiles: a and c identical, b opposite -> a, c adjacent
    m = pd.DataFrame({"a": [1, 2, 3, 4], "b": [4, 3, 2, 1], "c": [1, 2, 3, 4.1]})
    order = _cluster_order(m)
    assert abs(order.index("a") - order.index("c")) == 1
