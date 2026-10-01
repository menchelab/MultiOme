"""End-to-end on the shipped paper data (small subset) plus plot and CLI smoke tests."""

from __future__ import annotations

import matplotlib
import pytest

from multiome_algo import informed_rwr, layer_weights, modularity_table, significant_layers
from multiome_algo.cli import main
from multiome_algo.crossval import paper_configs, retrieval_cv
from multiome_algo.viz import (
    plot_candidate_ranking,
    plot_cv_performance,
    plot_layer_weights,
    plot_modularity_heatmap,
)
from multiome_core import filter_groups
from multiome_core.legacy import REPO_DATA, read_paper_groups, read_paper_multiplex

LAYERS = ["ppi", "HP", "coex_core"]
pytestmark = pytest.mark.skipif(not (REPO_DATA / "network_edgelists").exists(),
                                reason="paper data not available")


@pytest.fixture(scope="module")
def paper(tmp_path_factory):
    cache = tmp_path_factory.mktemp("cache")
    mpx = read_paper_multiplex(layer_ids=LAYERS, cache_dir=cache)
    groups = filter_groups(read_paper_groups())
    return mpx, groups


def test_paper_data_loads(paper):
    mpx, groups = paper
    assert mpx.layer_ids == LAYERS
    assert len(read_paper_groups()) == 28 and len(groups) == 26
    assert mpx["ppi"].adj.diagonal().sum() == 0
    assert "ppi" in mpx and mpx["HP"].tags


def test_pipeline_and_plots(paper, tmp_path):
    mpx, groups = paper
    gs = sorted(groups, key=len)[:3]
    tab = modularity_table(mpx, gs, n_trials=100)
    g = gs[0]
    zs = significant_layers(tab, g.id) or {"ppi": 1.0}
    seeds = sorted(g.genes)[: len(g.genes) // 2]
    res = informed_rwr(mpx, seeds, layer_weights=layer_weights(zs))
    assert res.table["rank"].notna().sum() > 1000
    cv = retrieval_cv(mpx, gs[:2], paper_configs(mpx), k_folds=3,
                      lcc_kwargs={"n_trials": 50})
    assert cv.folds["auroc"].dropna().between(0, 1).all()

    figs = [
        plot_modularity_heatmap(tab, layer_tags={lid: ";".join(mpx[lid].tags) for lid in LAYERS},
                                out=tmp_path / "h.png"),
        plot_layer_weights(tab, g.id, out=tmp_path / "w.png"),
        plot_candidate_ranking(res, highlight=set(g.genes), out=tmp_path / "r.png"),
        plot_cv_performance(cv, out=tmp_path / "cv.png"),
    ]
    assert all(isinstance(f, matplotlib.figure.Figure) for f in figs)
    assert len(list(tmp_path.glob("*.png"))) == 4


def test_cli_end_to_end(tmp_path, planted, capsys):
    mpx, groups = planted
    net = tmp_path / "net"
    net.mkdir()
    for layer in mpx:
        layer.edges()[["source", "target"]].to_csv(net / f"{layer.id}.tsv", sep="\t",
                                                   index=False)
    (tmp_path / "g.tsv").write_text(
        "group\tgene\n" + "".join(f"{g.id}\t{x}\n" for g in groups for x in sorted(g.genes))
    )
    common = ["--network", str(net), "--groups", str(tmp_path / "g.tsv"), "--no-cache",
              "--min-size", "10", "--trials", "50", "--out", str(tmp_path / "o")]
    main(["modularity", *common])
    main(["rank", *common, "--group", "module", "--top", "5"])
    main(["cv", *common, "--folds", "3", "--group", "module",
          "--configs", "informed=significant", "all=all", "c=C"])
    (tmp_path / "run.yaml").write_text(
        "command: modularity\n" + f"network: {net}\ngroups: {tmp_path / 'g.tsv'}\n"
        f"out: {tmp_path / 'o2'}\nmin_size: 10\ntrials: 20\nno_cache: true\nno_plots: true\n"
    )
    main(["run", str(tmp_path / "run.yaml")])
    out = tmp_path / "o"
    for f in ["modularity.tsv", "coverage.tsv", "modularity_heatmap.png", "ranking_module.tsv",
              "ranking_module.png", "cv_summary.tsv", "cv_folds.tsv", "cv_performance.png"]:
        assert (out / f).exists(), f
    assert (tmp_path / "o2" / "modularity.tsv").exists()
