"""Gene-group loaders. Each returns a list[GeneGroup] from some source."""

from multiome_algo.groups.orphanet import load_orphanet_groups

__all__ = ["load_orphanet_groups"]
