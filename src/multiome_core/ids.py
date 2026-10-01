"""Gene-symbol normalisation to current HGNC approved symbols.

Why: layers built at different times use different symbol versions. In the shipped 2021
multiplex the GTEx co-expression layers use post-2019 names (SEPTIN7, MARCHF5) while PPI,
GO, HPO, MPO and co-essentiality use the old ones (SEPT7, MARCH5), and the Orphanet gene
table has case variants (C9ORF72 vs C9orf72). Unnormalised, one gene becomes two nodes and
the multiplex coupling between its copies is lost.

Resolution order for each input symbol:
  1. exact current approved symbol
  2. unique previous symbol (``prev_symbol``)
  3. unique alias (``alias_symbol``)
  4. steps 1-3 case-insensitively
Ambiguous matches (a previous symbol/alias shared by several current genes) are left
unchanged and reported, never guessed.
"""

from __future__ import annotations

import urllib.request
from collections import defaultdict
from collections.abc import Iterable
from pathlib import Path

import pandas as pd

HGNC_URL = (
    "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
)
DEFAULT_PATH = Path.home() / ".cache" / "multiome" / "hgnc_complete_set.txt"


def load_hgnc(path: str | Path | None = None, download: bool = True) -> pd.DataFrame:
    """Load the HGNC complete set (columns symbol, prev_symbol, alias_symbol).

    Downloads to ``~/.cache/multiome/`` on first use unless `path` is given.
    """
    path = Path(path) if path else DEFAULT_PATH
    if not path.exists():
        if not download:
            raise FileNotFoundError(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        urllib.request.urlretrieve(HGNC_URL, path)  # noqa: S310
    df = pd.read_csv(path, sep="\t", dtype=str, usecols=["symbol", "prev_symbol", "alias_symbol"])
    return df.fillna("")


class SymbolMapper:
    """Map gene symbols to current HGNC approved symbols."""

    def __init__(self, hgnc: pd.DataFrame):
        self.current = set(hgnc["symbol"])
        self._prev: dict[str, set[str]] = defaultdict(set)
        self._alias: dict[str, set[str]] = defaultdict(set)
        for sym, prev, alias in hgnc[["symbol", "prev_symbol", "alias_symbol"]].itertuples(
            index=False
        ):
            for p in filter(None, prev.split("|")):
                self._prev[p].add(sym)
            for a in filter(None, alias.split("|")):
                self._alias[a].add(sym)
        self._upper: dict[str, set[str]] = defaultdict(set)
        for s in self.current:
            self._upper[s.upper()].add(s)
        self._prev_upper = self._fold(self._prev)
        self._alias_upper = self._fold(self._alias)

    @staticmethod
    def _fold(d: dict[str, set[str]]) -> dict[str, set[str]]:
        out: dict[str, set[str]] = defaultdict(set)
        for k, v in d.items():
            out[k.upper()] |= v
        return out

    @classmethod
    def default(cls, path: str | Path | None = None) -> SymbolMapper:
        return cls(load_hgnc(path))

    def resolve(self, symbol: str) -> tuple[str | None, str]:
        """Return (current_symbol or None, how) for one input symbol."""
        if symbol in self.current:
            return symbol, "current"
        u = symbol.upper()
        for table, how in (
            (self._prev, "previous"),
            (self._alias, "alias"),
            (self._upper, "case"),
            (self._prev_upper, "previous_case"),
            (self._alias_upper, "alias_case"),
        ):
            hits = table.get(u if how.endswith("case") else symbol)
            if hits:
                if len(hits) == 1:
                    return next(iter(hits)), how
                return None, "ambiguous:" + "|".join(sorted(hits))
        return None, "unmapped"

    def map(self, symbols: Iterable[str]) -> tuple[dict[str, str], pd.DataFrame]:
        """Map many symbols.

        Returns:
            mapping: {input: current} for every symbol that changes.
            report: one row per non-current input symbol (input, output, how).
        """
        mapping: dict[str, str] = {}
        rows = []
        for s in sorted(set(symbols)):
            out, how = self.resolve(s)
            if how == "current":
                continue
            if out is not None:
                mapping[s] = out
            rows.append({"input": s, "output": out or "", "how": how})
        return mapping, pd.DataFrame(rows, columns=["input", "output", "how"])


def normalize_symbols(multiplex, groups=(), mapper: SymbolMapper | None = None):
    """Normalise gene symbols in a Multiplex and gene groups to current HGNC symbols.

    Returns:
        (multiplex, groups, report): renamed copies plus the mapping report, which also
        lists unmapped/ambiguous symbols (left unchanged).
    """
    mapper = mapper or SymbolMapper.default()
    groups = list(groups)
    symbols = set(multiplex.universe)
    for g in groups:
        symbols |= g.genes
    mapping, report = mapper.map(symbols)
    return (
        multiplex.rename_genes(mapping) if mapping else multiplex,
        [g.rename(mapping) for g in groups],
        report,
    )
