"""Parity with the original R implementation (golden values from tests/r_parity/make_golden.R,
which runs the verbatim 2021 functions pmat_cal / supratransitional / RWR)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from multiome_algo.propagate import SupraOperator, informed_rwr, pmat_paper
from multiome_core.io import read_multiplex

HERE = Path(__file__).parent / "r_parity"
SEEDS = ["G02", "G07", "G11", "G27"]
CONFIGS = {"weighted": {"L1": 2.5, "L2": 4.0, "L3": 1.8}, "uniform": None}


@pytest.fixture(scope="module")
def toy():
    return read_multiplex(HERE / "toy", cache_dir=None)


@pytest.fixture(scope="module")
def golden():
    return pd.read_csv(HERE / "golden_rwr.tsv", sep="\t")


@pytest.mark.parametrize("config", CONFIGS)
def test_pmat_matches_r(config):
    w = list((CONFIGS[config] or dict.fromkeys(["L1", "L2", "L3"], 1.0)).values())
    expected = np.loadtxt(HERE / f"pmat_{config}.tsv")
    np.testing.assert_allclose(pmat_paper(w), expected, rtol=1e-12)


@pytest.mark.parametrize("config", CONFIGS)
def test_rwr_matches_r(toy, golden, config):
    exp = golden[golden["config"] == config].set_index("gene")
    res = informed_rwr(toy, SEEDS, layer_weights=CONFIGS[config], remove_seeds=False,
                       tol=1e-15)
    probs = res.layer_probs.loc[exp.index, ["L1", "L2", "L3"]]
    np.testing.assert_allclose(probs.to_numpy(), exp[["L1", "L2", "L3"]].to_numpy(),
                               rtol=1e-6, atol=1e-12)
    score = res.table.set_index("gene").loc[exp.index, "score"]
    np.testing.assert_allclose(score.to_numpy(), exp["avg"].to_numpy(), rtol=1e-6, atol=1e-12)


def test_operator_is_column_stochastic(toy):
    """Applying S to each unit vector gives columns summing to 1 (paper coupling) and the
    'present' variant conserves mass over present states."""
    for coupling in ("paper", "present"):
        op = SupraOperator(toy, toy.layer_ids, CONFIGS["weighted"], coupling=coupling)
        n, L = op.shape
        for k in range(n * L):
            g, layer = divmod(k, L)
            if coupling == "present" and not op.present[g, layer]:
                continue
            e = np.zeros((n, L))
            e[g, layer] = 1.0
            assert op.apply(e).sum() == pytest.approx(1.0)
