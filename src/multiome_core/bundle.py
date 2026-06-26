"""Read/write the multiplex bundle interchange format (.mpx).

A bundle is a directory:

    mybundle.mpx/
      manifest.yaml          # multiplex + per-layer metadata
      layers/
        <layer_id>.parquet   # edge table: source, target[, weight]

Parquet over 2-column TSV: typed, ~5-10x smaller, fast columnar reads on million-edge
layers, and carries weights + provenance without sidecar files. This is the only thing
the algorithm unit (multiome_algo) depends on; the network-generation unit (multiome_net)
produces it.
"""

from __future__ import annotations

from pathlib import Path

import networkx as nx
import pandas as pd
import yaml

from multiome_core.schema import BiologicalScale, Layer, Multiplex

MANIFEST = "manifest.yaml"
LAYERS_DIR = "layers"


def write_bundle(multiplex: Multiplex, path: str | Path) -> Path:
    """Write a Multiplex to a .mpx bundle directory. Returns the bundle path."""
    path = Path(path)
    layers_path = path / LAYERS_DIR
    layers_path.mkdir(parents=True, exist_ok=True)

    manifest: dict = {
        "name": multiplex.name,
        "node_namespace": multiplex.node_namespace,
        "layers": [],
    }
    for lid, layer in multiplex.layers.items():
        if layer.weighted:
            rows = [(u, v, d.get("weight", 1.0)) for u, v, d in layer.graph.edges(data=True)]
            df = pd.DataFrame(rows, columns=["source", "target", "weight"])
        else:
            df = pd.DataFrame(layer.graph.edges(), columns=["source", "target"])
        df.to_parquet(layers_path / f"{lid}.parquet", index=False)

        manifest["layers"].append(
            {
                "id": lid,
                "scale": str(layer.scale),
                "directed": layer.directed,
                "weighted": layer.weighted,
                "node_namespace": layer.node_namespace,
                "source": layer.source,
                "build_params": layer.build_params,
                "version": layer.version,
                "date": layer.date,
                "n_nodes": layer.n_nodes,
                "n_edges": layer.n_edges,
            }
        )

    (path / MANIFEST).write_text(yaml.safe_dump(manifest, sort_keys=False))
    return path


def read_bundle(path: str | Path, layer_ids: list[str] | None = None) -> Multiplex:
    """Read a .mpx bundle into a Multiplex. Optionally restrict to `layer_ids`."""
    path = Path(path)
    manifest = yaml.safe_load((path / MANIFEST).read_text())

    mpx = Multiplex(
        name=manifest.get("name", ""),
        node_namespace=manifest.get("node_namespace", "HGNC_symbol"),
    )
    for meta in manifest.get("layers", []):
        lid = meta["id"]
        if layer_ids is not None and lid not in layer_ids:
            continue
        df = pd.read_parquet(path / LAYERS_DIR / f"{lid}.parquet")
        directed = meta.get("directed", False)
        weighted = meta.get("weighted", False)
        create_using = nx.DiGraph if directed else nx.Graph
        if weighted and "weight" in df.columns:
            g = nx.from_pandas_edgelist(
                df, "source", "target", edge_attr="weight", create_using=create_using()
            )
        else:
            g = nx.from_pandas_edgelist(df, "source", "target", create_using=create_using())
        mpx.add(
            Layer(
                id=lid,
                graph=g,
                scale=BiologicalScale(meta.get("scale", "other")),
                directed=directed,
                weighted=weighted,
                node_namespace=meta.get("node_namespace", "HGNC_symbol"),
                source=meta.get("source", ""),
                build_params=meta.get("build_params", {}) or {},
                version=meta.get("version", ""),
                date=meta.get("date", ""),
            )
        )
    return mpx
