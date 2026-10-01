"""Shortcuts for the data shipped with the original 2021 repository.

``data/network_edgelists/`` holds the paper's 46 binary layers (header rows of varying
names - auto-detected), ``data/network_details.tsv`` their metadata, and
``data/table_disease_gene_assoc_orphanet_genetic.tsv`` the 28 Orphanet rare-genetic
disease groups (26 of which pass the paper's 20-2000 gene filter).
"""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

from multiome_core.io import DEFAULT_CACHE, read_gene_groups, read_multiplex
from multiome_core.schema import GeneGroup, Multiplex

REPO_DATA = Path(__file__).resolve().parents[2] / "data"


def read_paper_multiplex(
    data_dir: str | Path = REPO_DATA,
    layer_ids: Iterable[str] | None = None,
    cache_dir: str | Path | None = DEFAULT_CACHE,
) -> Multiplex:
    """Load the paper's multiplex (all 46 layers unless `layer_ids` is given)."""
    data_dir = Path(data_dir)
    return read_multiplex(
        data_dir / "network_edgelists",
        layer_ids=layer_ids,
        metadata=data_dir / "network_details.tsv",
        cache_dir=cache_dir,
        name="buphamalai2021",
    )


def read_paper_groups(data_dir: str | Path = REPO_DATA) -> list[GeneGroup]:
    """Load the 28 Orphanet rare-genetic-disease groups (unfiltered)."""
    path = Path(data_dir) / "table_disease_gene_assoc_orphanet_genetic.tsv"
    return [
        GeneGroup(id=g.id, genes=g.genes, label=g.label, source="orphanet")
        for g in read_gene_groups(path)
    ]


def read_paper_dataset(
    data_dir: str | Path = REPO_DATA,
    layer_ids: Iterable[str] | None = None,
    normalize: bool = True,
    cache_dir: str | Path | None = DEFAULT_CACHE,
):
    """The paper's multiplex and gene groups, by default with symbols normalised to
    current HGNC (the shipped layers mix old and new symbol versions; see DEVIATIONS).

    Returns:
        (multiplex, groups, report): `report` lists every renamed/unresolved symbol
        (empty DataFrame when `normalize` is False).
    """
    import pandas as pd

    mpx = read_paper_multiplex(data_dir, layer_ids, cache_dir)
    groups = read_paper_groups(data_dir)
    if not normalize:
        return mpx, groups, pd.DataFrame(columns=["input", "output", "how"])
    from multiome_core.ids import normalize_symbols

    return normalize_symbols(mpx, groups)
