from __future__ import annotations

import numpy as np

from multiome_algo.crossval import CVConfig, make_folds, retrieval_cv
from multiome_core.schema import GeneGroup


def test_make_folds_deterministic_and_partitioning():
    genes = [f"g{i}" for i in range(23)]
    f1 = make_folds(genes, 5, seed=0, group_id="x")
    assert f1 == make_folds(reversed(genes), 5, seed=0, group_id="x")
    assert sorted(sum(f1, [])) == sorted(genes)
    assert {len(f) for f in f1} <= {4, 5}
    assert f1 != make_folds(genes, 5, seed=1, group_id="x")


def _configs():
    return [CVConfig("informed", "significant", "z"), CVConfig("all_uniform", "all", "uniform"),
            CVConfig("C_only", ["C"], "uniform")]


def test_cv_on_planted_module(planted):
    mpx, groups = planted
    for protocol in ("train", "paper"):
        cv = retrieval_cv(mpx, groups[:1], _configs(), k_folds=5, protocol=protocol,
                          min_size=None, lcc_kwargs={"n_trials": 100})
        s = cv.summary().set_index("config")
        assert s.loc["informed", "auroc_median"] > 0.9
        assert s.loc["informed", "auroc_median"] > s.loc["C_only", "auroc_median"]
        assert (s["n_folds"] == 5).all()
        # folds are paired across configs
        held = cv.held_out.groupby(["config", "fold"])["gene"].apply(frozenset).unstack(0)
        assert (held["informed"] == held["all_uniform"]).all()
        assert set(cv.roc["fpr"].round(6)) == set(np.linspace(0, 1, 101).round(6))
        n_tables = 5 if protocol == "train" else 1
        assert len(cv.modularity) == n_tables


def test_cv_no_layers_status(planted):
    mpx, _ = planted
    # 11 genes, 5 folds -> <= 9 training genes per layer, below min_genes=10: never assessed
    small = GeneGroup("small", [f"G{i:03d}" for i in range(11)])
    cv = retrieval_cv(mpx, [small], _configs()[:1], k_folds=5, min_size=None,
                      lcc_kwargs={"n_trials": 50})
    assert (cv.folds["status"] == "no_layers").all()
    assert cv.folds["auroc"].isna().all()
