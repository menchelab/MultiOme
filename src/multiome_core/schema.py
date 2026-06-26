"""Typed data model for the multiplex network interchange contract.

A `Multiplex` is a set of `Layer`s over a shared node namespace (genes). A `GeneGroup`
is any set of genes (a rare-disease gene set, a PanelApp panel, a GO term's genes, or a
custom list) — the generic unit of analysis. The algorithm (multiome_algo) operates on
these objects and is agnostic to how the networks were built.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import StrEnum

import networkx as nx


class BiologicalScale(StrEnum):
    """The level of biological organization a layer represents (Buphamalai et al. 2021)."""

    GENOME = "genome"  # e.g. CRISPR co-essentiality
    TRANSCRIPTOME = "transcriptome"  # co-expression (tissue-specific + core)
    PROTEOME = "proteome"  # protein-protein interactions
    PATHWAY = "pathway"  # pathway co-membership
    FUNCTION = "function"  # GO biological-process / molecular-function similarity
    PHENOTYPE = "phenotype"  # HPO / MPO phenotypic similarity
    OTHER = "other"


@dataclass
class Layer:
    """A single network layer over the shared gene namespace.

    Attributes:
        id: Unique layer identifier (e.g. "ppi", "coex_BRC").
        graph: The network. Nodes are gene identifiers in `node_namespace`.
        scale: Biological scale this layer represents.
        directed: Whether edges are directed.
        weighted: Whether edges carry meaningful weights (stored on the "weight" edge attr).
        node_namespace: Identifier system for nodes (e.g. "HGNC_symbol").
        source: Free-text provenance (data source / paper / file).
        build_params: How the layer was constructed (thresholds, method, etc.).
        version: Source/build version string.
        date: Build or retrieval date (ISO string).
    """

    id: str
    graph: nx.Graph
    scale: BiologicalScale = BiologicalScale.OTHER
    directed: bool = False
    weighted: bool = False
    node_namespace: str = "HGNC_symbol"
    source: str = ""
    build_params: dict = field(default_factory=dict)
    version: str = ""
    date: str = ""

    @property
    def n_nodes(self) -> int:
        return self.graph.number_of_nodes()

    @property
    def n_edges(self) -> int:
        return self.graph.number_of_edges()

    def nodes(self) -> set[str]:
        return set(self.graph.nodes())


@dataclass
class Multiplex:
    """A collection of layers over a shared gene namespace, with metadata.

    Layers need not span identical node sets; the universe is the union of all layer nodes.
    """

    layers: dict[str, Layer] = field(default_factory=dict)
    name: str = ""
    node_namespace: str = "HGNC_symbol"

    def add(self, layer: Layer) -> None:
        if layer.id in self.layers:
            raise ValueError(f"Duplicate layer id: {layer.id!r}")
        self.layers[layer.id] = layer

    @property
    def layer_ids(self) -> list[str]:
        return list(self.layers.keys())

    def universe(self) -> set[str]:
        """Union of all genes appearing in any layer."""
        u: set[str] = set()
        for layer in self.layers.values():
            u |= layer.nodes()
        return u

    def __len__(self) -> int:
        return len(self.layers)

    def __getitem__(self, layer_id: str) -> Layer:
        return self.layers[layer_id]


@dataclass
class GeneGroup:
    """A set of genes — the generic unit of analysis (need not be a disease).

    Attributes:
        id: Stable identifier.
        label: Human-readable name.
        genes: The gene set, in the same namespace as the multiplex.
        source: Provenance (e.g. "orphanet", "panelapp", "go", "custom").
        metadata: Anything extra (ontology id, inheritance, etc.).
    """

    id: str
    genes: set[str]
    label: str = ""
    source: str = ""
    metadata: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        # Accept any iterable of genes; store as a set.
        if not isinstance(self.genes, set):
            self.genes = set(self.genes)

    def __len__(self) -> int:
        return len(self.genes)

    def in_universe(self, universe: set[str]) -> set[str]:
        """Genes of this group that are present in the given node universe."""
        return self.genes & universe
