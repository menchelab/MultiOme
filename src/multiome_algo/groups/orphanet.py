"""Load the legacy OrphaNet disease-gene-association table into GeneGroups.

Format: TSV with columns `name` and `all_genes`, where all_genes is a
semicolon-separated list of HGNC symbols. Each row becomes one GeneGroup.
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

from multiome_core.schema import GeneGroup


def _slugify(name: str) -> str:
    s = re.sub(r"[^0-9a-zA-Z]+", "_", name.strip().lower())
    return s.strip("_")


def load_orphanet_groups(path: str | Path) -> list[GeneGroup]:
    """Read the disease-gene-association table into a list of GeneGroups."""
    path = Path(path)
    df = pd.read_csv(path, sep="\t", dtype=str).dropna(subset=["name", "all_genes"])

    groups: list[GeneGroup] = []
    for _, row in df.iterrows():
        genes = {g.strip() for g in row["all_genes"].split(";") if g.strip()}
        groups.append(
            GeneGroup(
                id=_slugify(row["name"]),
                genes=genes,
                label=row["name"],
                source="orphanet",
            )
        )
    return groups
