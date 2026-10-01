"""Typed data model: network layers, the multiplex, and gene groups.

A `Layer` stores one undirected network as a symmetric sparse adjacency matrix over its
own gene index. A `Multiplex` is a set of layers over a shared gene namespace; it can
align every layer to a common universe index, which is what propagation operates on.
A `GeneGroup` is any gene set (a rare-disease gene set, a clinical panel, a GO term, a
custom list) - the generic unit of analysis.

Sparse matrices rather than networkx graphs: the shipped 46-layer multiplex has ~21M
edges, which networkx cannot hold or traverse at a usable speed. `Layer.to_networkx()`
is available for inspection.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass, field
from functools import cached_property

import numpy as np
import pandas as pd
import scipy.sparse as sp


@dataclass(eq=False)
class Layer:
    """One undirected network layer.

    Attributes:
        id: Unique layer identifier (e.g. "ppi", "coex_BRC").
        genes: Gene identifiers, one per row/column of `adj` (sorted, unique).
        adj: Symmetric CSR adjacency (n x n). Binary unless `weighted`; no self-loops.
        weighted: Whether `adj` carries meaningful edge weights.
        tags: Free-text labels for grouping layers in plots/tables (e.g. "co-expression").
        description: Human-readable description (e.g. tissue name).
        source: Provenance (file path, database, ...).
    """

    id: str
    genes: np.ndarray
    adj: sp.csr_matrix
    weighted: bool = False
    tags: tuple[str, ...] = ()
    description: str = ""
    source: str = ""

    def __post_init__(self) -> None:
        self.genes = np.asarray(self.genes, dtype=object)
        if self.adj.shape != (len(self.genes), len(self.genes)):
            raise ValueError(
                f"Layer {self.id!r}: adjacency shape {self.adj.shape} does not match "
                f"{len(self.genes)} genes."
            )
        self.adj = sp.csr_matrix(self.adj)

    @classmethod
    def from_edges(
        cls,
        layer_id: str,
        sources: Iterable[str],
        targets: Iterable[str],
        weights: Iterable[float] | None = None,
        **meta,
    ) -> Layer:
        """Build a layer from parallel source/target (and optional weight) sequences.

        Self-loops and missing values are dropped; duplicate (undirected) edges are
        collapsed, keeping the maximum weight.
        """
        df = pd.DataFrame({"source": list(sources), "target": list(targets)})
        df["weight"] = 1.0 if weights is None else np.asarray(list(weights), dtype=float)
        df = df.dropna()
        df["source"] = df["source"].astype(str)
        df["target"] = df["target"].astype(str)
        df = df[df["source"] != df["target"]]

        genes = np.array(sorted(set(df["source"]) | set(df["target"])), dtype=object)
        idx = pd.Index(genes)
        i = idx.get_indexer(df["source"])
        j = idx.get_indexer(df["target"])
        lo, hi = np.minimum(i, j), np.maximum(i, j)
        edges = (
            pd.DataFrame({"lo": lo, "hi": hi, "w": df["weight"].to_numpy()})
            .groupby(["lo", "hi"], sort=False)["w"]
            .max()
            .reset_index()
        )
        n = len(genes)
        rows = np.concatenate([edges["lo"], edges["hi"]])
        cols = np.concatenate([edges["hi"], edges["lo"]])
        vals = np.concatenate([edges["w"], edges["w"]])
        adj = sp.csr_matrix((vals, (rows, cols)), shape=(n, n))
        return cls(id=layer_id, genes=genes, adj=adj, weighted=weights is not None, **meta)

    @property
    def n_nodes(self) -> int:
        return len(self.genes)

    @property
    def n_edges(self) -> int:
        return self.adj.nnz // 2

    @cached_property
    def index(self) -> pd.Index:
        """Gene -> row position."""
        return pd.Index(self.genes)

    @cached_property
    def gene_set(self) -> frozenset[str]:
        return frozenset(self.genes)

    @cached_property
    def degree(self) -> np.ndarray:
        """Unweighted node degree."""
        return np.diff(self.adj.indptr)

    def positions(self, genes: Iterable[str]) -> np.ndarray:
        """Row positions of the given genes that are present in this layer."""
        pos = self.index.get_indexer(list(genes))
        return pos[pos >= 0]

    def edges(self) -> pd.DataFrame:
        """Edge table (each undirected edge once)."""
        upper = sp.triu(self.adj, k=1).tocoo()
        return pd.DataFrame(
            {
                "source": self.genes[upper.row],
                "target": self.genes[upper.col],
                "weight": upper.data,
            }
        )

    def rename(self, mapping: dict[str, str]) -> Layer:
        """Return a copy with genes renamed; genes that collapse onto one name are merged.
        Returns `self` when no gene of this layer is renamed."""
        if not any(g in mapping for g in self.genes):
            return self
        e = self.edges()
        return Layer.from_edges(
            self.id,
            e["source"].map(lambda g: mapping.get(g, g)),
            e["target"].map(lambda g: mapping.get(g, g)),
            e["weight"] if self.weighted else None,
            tags=self.tags,
            description=self.description,
            source=self.source,
        )

    def to_networkx(self):
        import networkx as nx

        e = self.edges()
        attr = "weight" if self.weighted else None
        return nx.from_pandas_edgelist(e, "source", "target", edge_attr=attr)


@dataclass(eq=False)
class Multiplex:
    """A collection of layers over a shared gene namespace.

    Layers need not span identical gene sets; the universe is the union of all layer genes.
    """

    layers: dict[str, Layer] = field(default_factory=dict)
    name: str = ""

    def add(self, layer: Layer) -> None:
        if layer.id in self.layers:
            raise ValueError(f"Duplicate layer id: {layer.id!r}")
        self.layers[layer.id] = layer
        self._invalidate()

    def _invalidate(self) -> None:
        self.__dict__.pop("universe", None)
        self.__dict__.pop("universe_index", None)
        self.__dict__.pop("_aligned", None)
        self.__dict__.pop("_normalized", None)

    @property
    def layer_ids(self) -> list[str]:
        return list(self.layers)

    @cached_property
    def universe(self) -> np.ndarray:
        """Sorted union of all genes appearing in any layer."""
        genes: set[str] = set()
        for layer in self.layers.values():
            genes |= layer.gene_set
        return np.array(sorted(genes), dtype=object)

    @cached_property
    def universe_index(self) -> pd.Index:
        return pd.Index(self.universe)

    def gene_set(self) -> frozenset[str]:
        return frozenset(self.universe)

    @cached_property
    def _aligned(self) -> dict[str, sp.csr_matrix]:
        return {}

    def aligned_adjacency(self, layer_id: str) -> sp.csr_matrix:
        """Layer adjacency re-indexed onto the multiplex universe (N x N), cached."""
        if layer_id not in self._aligned:
            layer = self.layers[layer_id]
            pos = self.universe_index.get_indexer(layer.genes)
            coo = layer.adj.tocoo()
            n = len(self.universe)
            self._aligned[layer_id] = sp.csr_matrix(
                (coo.data, (pos[coo.row], pos[coo.col])), shape=(n, n)
            )
        return self._aligned[layer_id]

    def presence(self, layer_ids: list[str] | None = None) -> np.ndarray:
        """Boolean matrix (N genes x L layers): is gene present in layer."""
        layer_ids = layer_ids or self.layer_ids
        mask = np.zeros((len(self.universe), len(layer_ids)), dtype=bool)
        for k, lid in enumerate(layer_ids):
            mask[self.universe_index.get_indexer(self.layers[lid].genes), k] = True
        return mask

    def subset(self, layer_ids: Iterable[str]) -> Multiplex:
        """A new Multiplex over the given layers (layers are shared, not copied)."""
        return Multiplex({lid: self.layers[lid] for lid in layer_ids}, name=self.name)

    def rename_genes(self, mapping: dict[str, str]) -> Multiplex:
        """Apply a gene-id mapping to every layer (see multiome_core.ids)."""
        return Multiplex(
            {lid: layer.rename(mapping) for lid, layer in self.layers.items()}, name=self.name
        )

    def summary(self) -> pd.DataFrame:
        return pd.DataFrame(
            [
                {
                    "layer_id": lid,
                    "n_nodes": layer.n_nodes,
                    "n_edges": layer.n_edges,
                    "weighted": layer.weighted,
                    "tags": ";".join(layer.tags),
                    "description": layer.description,
                }
                for lid, layer in self.layers.items()
            ]
        )

    def __len__(self) -> int:
        return len(self.layers)

    def __getitem__(self, layer_id: str) -> Layer:
        return self.layers[layer_id]

    def __iter__(self):
        return iter(self.layers.values())

    def __contains__(self, layer_id: str) -> bool:
        return layer_id in self.layers


@dataclass(eq=False)
class GeneGroup:
    """A set of genes - the generic unit of analysis (need not be a disease).

    Attributes:
        id: Stable identifier (taken from the source where possible, e.g. an ontology id).
        genes: The gene set, in the same namespace as the multiplex.
        label: Human-readable name.
        source: Provenance (e.g. "orphanet", "gmt:<file>", "custom").
        metadata: Anything extra.
    """

    id: str
    genes: frozenset[str]
    label: str = ""
    source: str = ""
    metadata: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.genes = frozenset(self.genes)
        self.label = self.label or self.id

    def __len__(self) -> int:
        return len(self.genes)

    def in_universe(self, universe: Iterable[str]) -> frozenset[str]:
        """Genes of this group that are present in the given universe."""
        return self.genes & frozenset(universe)

    def rename(self, mapping: dict[str, str]) -> GeneGroup:
        return GeneGroup(
            id=self.id,
            genes=frozenset(mapping.get(g, g) for g in self.genes),
            label=self.label,
            source=self.source,
            metadata=dict(self.metadata),
        )


def filter_groups(
    groups: Iterable[GeneGroup],
    min_size: int | None = 20,
    max_size: int | None = 2000,
    universe: Iterable[str] | None = None,
) -> list[GeneGroup]:
    """Keep groups whose size lies in [min_size, max_size].

    Size is the total number of annotated genes (the paper's rule, 20-2000) unless a
    `universe` is given, in which case only genes present in it are counted.
    """
    uni = frozenset(universe) if universe is not None else None
    out = []
    for g in groups:
        n = len(g.genes & uni) if uni is not None else len(g)
        if (min_size is None or n >= min_size) and (max_size is None or n <= max_size):
            out.append(g)
    return out
