from __future__ import annotations

import pandas as pd

from multiome_core.ids import SymbolMapper, normalize_symbols
from multiome_core.schema import GeneGroup, Layer, Multiplex

HGNC = pd.DataFrame(
    {
        "symbol": ["SEPTIN7", "MARCHF5", "C9orf72", "GENEA", "GENEB"],
        "prev_symbol": ["SEPT7", "MARCH5", "", "SHARED", "SHARED"],
        "alias_symbol": ["CDC10", "", "ALS-FTD", "", ""],
    }
)


def test_mapper_resolution_order():
    m = SymbolMapper(HGNC)
    mapping, report = m.map(["SEPTIN7", "SEPT7", "CDC10", "C9ORF72", "SHARED", "XYZ", "sept7"])
    assert mapping == {"SEPT7": "SEPTIN7", "CDC10": "SEPTIN7", "C9ORF72": "C9orf72",
                       "sept7": "SEPTIN7"}
    how = dict(zip(report["input"], report["how"], strict=True))
    assert how["SEPT7"] == "previous" and how["CDC10"] == "alias"
    assert how["C9ORF72"] == "case"
    assert how["SHARED"].startswith("ambiguous")  # never guessed
    assert how["XYZ"] == "unmapped"
    assert "SEPTIN7" not in how  # current symbols are not reported


def test_normalize_merges_versions():
    # the same gene under old and new symbols in two layers becomes one node
    mpx = Multiplex()
    mpx.add(Layer.from_edges("old", ["SEPT7"], ["MARCH5"]))
    mpx.add(Layer.from_edges("new", ["SEPTIN7"], ["GENEA"]))
    assert len(mpx.universe) == 4
    out, groups, report = normalize_symbols(mpx, [GeneGroup("g", ["SEPT7", "GENEA"])],
                                            mapper=SymbolMapper(HGNC))
    assert set(out.universe) == {"SEPTIN7", "MARCHF5", "GENEA"}
    assert groups[0].genes == {"SEPTIN7", "GENEA"}
    assert out["new"] is mpx["new"]  # untouched layers are reused
    assert len(report) == 2
