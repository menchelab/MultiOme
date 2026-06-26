"""Read the original MultiOme TSV edge lists into the typed Layer/Multiplex model.

The 2021 repo ships 46 two-column TSV edge lists in data/network_edgelists/ with
inconsistent headers ("A"/"B", "name1"/"name2", "Gene Name Interactor A/B", ...).
We ignore the header names and take the first two columns as (source, target).

This is the phase-1 bridge: it lets the algorithm (multiome_algo) run against the
paper's own networks before the network-generation unit (multiome_net) is built.
"""

from __future__ import annotations

from pathlib import Path

import networkx as nx
import pandas as pd

from multiome_core.schema import BiologicalScale, Layer, Multiplex

# Layer id -> biological scale for the original 46-layer multiplex.
# Co-expression layers (coex_*) are mapped programmatically below.
_SCALE_BY_ID: dict[str, BiologicalScale] = {
    "ppi": BiologicalScale.PROTEOME,
    "co-essential": BiologicalScale.GENOME,
    "reactome_copathway": BiologicalScale.PATHWAY,
    "GOBP": BiologicalScale.FUNCTION,
    "GOMF": BiologicalScale.FUNCTION,
    "HP": BiologicalScale.PHENOTYPE,
    "MP": BiologicalScale.PHENOTYPE,
}


def _scale_for(layer_id: str) -> BiologicalScale:
    if layer_id.startswith("coex_"):
        return BiologicalScale.TRANSCRIPTOME
    return _SCALE_BY_ID.get(layer_id, BiologicalScale.OTHER)


def read_edgelist(path: str | Path, layer_id: str | None = None) -> Layer:
    """Read one TSV edge list into a Layer.

    The first two columns are treated as (source, target) regardless of header names.
    Self-loops are dropped. Graph is undirected and unweighted (matches the originals).
    """
    path = Path(path)
    layer_id = layer_id or path.stem
    df = pd.read_csv(path, sep="\t", usecols=[0, 1], header=0, dtype=str)
    df.columns = ["source", "target"]
    df = df.dropna()
    df = df[df["source"] != df["target"]]  # drop self-loops

    g = nx.from_pandas_edgelist(df, "source", "target")
    return Layer(
        id=layer_id,
        graph=g,
        scale=_scale_for(layer_id),
        directed=False,
        weighted=False,
        node_namespace="HGNC_symbol",
        source=f"legacy:{path.name}",
    )


def read_multiplex(
    edgelist_dir: str | Path,
    layer_ids: list[str] | None = None,
    name: str = "legacy",
) -> Multiplex:
    """Read a directory of TSV edge lists into a Multiplex.

    Args:
        edgelist_dir: Directory containing *.tsv edge lists.
        layer_ids: If given, only load these layer ids (file stems). Otherwise load all *.tsv.
        name: Name for the multiplex.
    """
    edgelist_dir = Path(edgelist_dir)
    if layer_ids is None:
        paths = sorted(edgelist_dir.glob("*.tsv"))
    else:
        paths = [edgelist_dir / f"{lid}.tsv" for lid in layer_ids]
        missing = [p for p in paths if not p.exists()]
        if missing:
            raise FileNotFoundError(f"Missing edge lists: {[str(p) for p in missing]}")

    mpx = Multiplex(name=name)
    for p in paths:
        mpx.add(read_edgelist(p))
    return mpx
