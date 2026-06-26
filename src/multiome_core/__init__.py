"""multiome_core — the interchange contract shared by network generation and the algorithm.

Defines the typed data model (Layer, Multiplex, GeneGroup), the on-disk multiplex
bundle format (.mpx), and a reader for the legacy TSV edge lists. Both multiome_net
(produces bundles) and multiome_algo (consumes bundles) depend only on this package.
"""

from multiome_core.schema import BiologicalScale, GeneGroup, Layer, Multiplex

__all__ = ["BiologicalScale", "GeneGroup", "Layer", "Multiplex"]
