from __future__ import annotations

import numpy as np
import pytest

from multiome_core import coverage_report, read_edgelist, read_gene_groups, read_multiplex
from multiome_core.io import write_gmt
from multiome_core.schema import GeneGroup, Layer


def test_from_edges_cleans(tmp_path):
    layer = Layer.from_edges("x", ["A", "B", "A", "C", "A", None], ["B", "A", "A", "D", "C", "D"])
    # self-loop A-A, NA and duplicate B-A dropped
    assert layer.n_edges == 3
    assert (layer.adj != layer.adj.T).nnz == 0
    assert layer.adj.diagonal().sum() == 0
    assert list(layer.genes) == ["A", "B", "C", "D"]


@pytest.mark.parametrize(
    "content,sep_name",
    [
        ("gene1\tgene2\nA\tB\nB\tC\n", "tsv-header"),
        ("A\tB\nB\tC\nC\tD\n", "tsv-noheader"),
        ("Gene Name A,Gene Name B\nA,B\nB,C\n", "csv-spaced-header"),
        ('A\t"NA"\nA\tB\nB\tC\n', "quoted-na"),
    ],
)
def test_read_edgelist_header_and_delimiter(tmp_path, content, sep_name):
    f = tmp_path / "L.tsv"
    f.write_text(content)
    layer = read_edgelist(f)
    assert layer.id == "L"
    assert set(layer.genes) <= {"A", "B", "C", "D"}
    assert layer.n_edges >= 2


def test_read_edgelist_weights(tmp_path):
    f = tmp_path / "w.csv"
    f.write_text("a,b,w\nA,B,0.5\nB,A,2.0\nB,C,1\n")
    layer = read_edgelist(f, weights=True)
    assert layer.weighted
    i, j = layer.index.get_indexer(["A", "B"])
    assert layer.adj[i, j] == 2.0  # duplicate undirected edge keeps the max


def test_read_multiplex_metadata_and_cache(tmp_path):
    net = tmp_path / "net"
    net.mkdir()
    (net / "ppi.tsv").write_text("A\tB\nB\tC\n")
    (net / "coex.csv").write_text("A,C\nC,D\n")
    (net / "layers.tsv").write_text("layer_id\ttags\tdescription\nppi\tmolecular\tPPI\n")
    cache = tmp_path / "cache"
    m1 = read_multiplex(net, cache_dir=cache)
    m2 = read_multiplex(net, cache_dir=cache)  # from cache
    assert m1.layer_ids == ["coex", "ppi"]
    assert m1["ppi"].tags == ("molecular",) and m1["ppi"].description == "PPI"
    assert list(m1.universe) == ["A", "B", "C", "D"]
    for lid in m1.layer_ids:
        assert (m1[lid].adj != m2[lid].adj).nnz == 0
        assert list(m1[lid].genes) == list(m2[lid].genes)
    assert len(list(cache.iterdir())) == 2
    with pytest.raises(FileNotFoundError):
        read_multiplex(net, layer_ids=["nope"], cache_dir=None)


def test_gene_group_formats(tmp_path):
    gmt = tmp_path / "s.gmt"
    write_gmt([GeneGroup("S1", ["A", "B"], label="set one")], gmt)
    (g,) = read_gene_groups(gmt)
    assert g.id == "S1" and g.genes == {"A", "B"}

    long = tmp_path / "long.tsv"
    long.write_text("group\tgene\nS1\tA\nS1\tB\nS2\tC\n")
    gs = {g.id: g.genes for g in read_gene_groups(long)}
    assert gs == {"S1": {"A", "B"}, "S2": {"C"}}

    wide = tmp_path / "wide.tsv"
    wide.write_text("orphaID\tname\tall_genes\n1\tDis\tA;B, C\n2\tDis\tD|E\n2\tX\tF\n")
    gs = read_gene_groups(wide)
    assert [g.id for g in gs] == ["1", "2", "2#1"]
    assert gs[0].label == "Dis" and gs[0].genes == {"A", "B", "C"}
    assert gs[1].genes == {"D", "E"}


def test_coverage_report(planted):
    mpx, groups = planted
    g = GeneGroup("x", ["G000", "G260", "NOPE"])
    cov = coverage_report(mpx, [g]).iloc[0]
    assert cov["n_genes"] == 3 and cov["n_in_universe"] == 2 and cov["n_missing"] == 1
    assert cov["n_in:A"] == 2 and cov["n_in:B"] == 1
    assert "NOPE" in cov["missing"]


def test_presence_and_alignment(planted):
    mpx, _ = planted
    pres = mpx.presence()
    assert pres.shape == (len(mpx.universe), 3)
    assert not pres[mpx.universe_index.get_loc("G260"), 1]
    a = mpx.aligned_adjacency("A")
    assert a.shape == (len(mpx.universe),) * 2
    assert a.nnz == mpx["A"].adj.nnz
    assert np.isclose(a.sum(), mpx["A"].adj.sum())
