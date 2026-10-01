"""multiome_core - data model and input handling.

Typed model (Layer, Multiplex, GeneGroup), readers for edge-list folders and gene-set
files, HGNC symbol normalisation, and shortcuts for the data shipped with the paper.
"""

from multiome_core.io import (
    coverage_report,
    multiplex_from_dataframes,
    read_edgelist,
    read_gene_groups,
    read_multiplex,
    write_gmt,
)
from multiome_core.schema import GeneGroup, Layer, Multiplex, filter_groups

__all__ = [
    "GeneGroup",
    "Layer",
    "Multiplex",
    "coverage_report",
    "filter_groups",
    "multiplex_from_dataframes",
    "read_edgelist",
    "read_gene_groups",
    "read_multiplex",
    "write_gmt",
]
