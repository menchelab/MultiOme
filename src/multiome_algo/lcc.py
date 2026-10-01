"""Largest-connected-component (LCC) modularity of a gene group in each layer.

Port of ``LCC_functions.R`` / ``process_LCC_result.R`` (Buphamalai et al. 2021). For a
group with n genes present in a layer, the observed LCC size of the induced subgraph is
compared with the LCC sizes of random gene sets of the same size drawn from the layer:

    z = (LCC_obs - mean_rand) / sd_rand          (sd with ddof=1, as R's sd())
    p_norm = 1 - Phi(z)                          (paper)
    p_emp  = (#{rand >= obs} + 1) / (trials + 1) (permutation p; option)

`modularity_table` applies the paper's significance rule: BH correction over *all*
group x layer pairs in the table, significant if q < 0.05 with >= 10 group genes in the
layer and an LCC of >= 5.

Nulls
    "uniform" (paper): random sets drawn uniformly from the layer's genes.
    "degree": degree-preserving - each group gene is replaced by a random gene from the
        same degree bin. Corrects for well-studied (high-degree) genes forming large LCCs
        by chance, which inflates z under the uniform null.

Null samples are cached per (layer, set size) - or per degree composition - and drawn
from an RNG seeded by (seed, layer id, size), so results do not depend on the order in
which groups/layers are processed.
"""

from __future__ import annotations

import zlib
from collections.abc import Iterable
from dataclasses import asdict, dataclass

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components
from scipy.stats import norm

from multiome_core.schema import GeneGroup, Layer, Multiplex

NULLS = ("uniform", "degree")


@dataclass
class LCCResult:
    group_id: str
    layer_id: str
    n_genes_in_layer: int
    lcc_size: int
    rand_mean: float
    rand_std: float
    z_score: float
    p_norm: float
    p_emp: float


def lcc_size(adj: sp.csr_matrix, positions: np.ndarray) -> int:
    """Size of the largest connected component of the subgraph induced by `positions`."""
    if len(positions) == 0:
        return 0
    sub = adj[positions][:, positions]
    _, labels = connected_components(sub, directed=False)
    return int(np.bincount(labels).max())


def _stable_hash(*parts) -> int:
    return zlib.crc32("|".join(map(str, parts)).encode())


def degree_bins(degree: np.ndarray, min_bin_size: int = 100) -> np.ndarray:
    """Assign each node a degree bin; nodes of equal degree share a bin, and every bin
    holds at least `min_bin_size` nodes (the last bin absorbs any remainder)."""
    order = np.argsort(degree, kind="stable")
    bins = np.empty(len(degree), dtype=int)
    b, count, i = 0, 0, 0
    n = len(degree)
    while i < n:
        j = i
        while j < n and degree[order[j]] == degree[order[i]]:
            j += 1
        bins[order[i:j]] = b
        count += j - i
        i = j
        if count >= min_bin_size and n - i >= min_bin_size:
            b += 1
            count = 0
    return bins


class LCCNull:
    """Random-set LCC-size distributions for one layer, cached by set size/composition."""

    def __init__(
        self,
        layer: Layer,
        n_trials: int = 1000,
        null: str = "uniform",
        seed: int = 0,
        min_bin_size: int = 100,
    ):
        if null not in NULLS:
            raise ValueError(f"null must be one of {NULLS}, got {null!r}")
        self.layer = layer
        self.n_trials = n_trials
        self.null = null
        self.seed = seed
        self._cache: dict[tuple, np.ndarray] = {}
        if null == "degree":
            self.bins = degree_bins(layer.degree, min_bin_size)
            self._members = [np.flatnonzero(self.bins == b) for b in range(self.bins.max() + 1)]

    def samples(self, positions: np.ndarray) -> np.ndarray:
        """LCC sizes of `n_trials` random sets matched to the given group positions."""
        if self.null == "uniform":
            key = (len(positions),)
        else:
            key = tuple(np.bincount(self.bins[positions], minlength=len(self._members)))
        if key not in self._cache:
            self._cache[key] = self._draw(key)
        return self._cache[key]

    def _draw(self, key: tuple) -> np.ndarray:
        rng = np.random.default_rng([self.seed, _stable_hash(self.layer.id, self.null, *key)])
        adj = self.layer.adj
        out = np.empty(self.n_trials)
        if self.null == "uniform":
            (n,) = key
            for t in range(self.n_trials):
                out[t] = lcc_size(adj, rng.choice(self.layer.n_nodes, size=n, replace=False))
        else:
            for t in range(self.n_trials):
                pos = np.concatenate(
                    [rng.choice(m, size=c, replace=False)
                     for m, c in zip(self._members, key, strict=True) if c]
                )
                out[t] = lcc_size(adj, pos)
        return out


def lcc_modularity(
    layer: Layer,
    group: GeneGroup | Iterable[str],
    n_trials: int = 1000,
    min_genes: int = 10,
    null: str | LCCNull = "uniform",
    seed: int = 0,
    group_id: str | None = None,
) -> LCCResult | None:
    """LCC z-score of a gene group in one layer.

    Returns None if fewer than `min_genes` group genes are in the layer (the paper's
    ``minnode = 10``). Pass a shared `LCCNull` to reuse null samples across groups.
    """
    if isinstance(group, GeneGroup):
        group_id = group_id or group.id
        genes = group.genes
    else:
        genes = frozenset(group)
    pos = layer.positions(genes)
    if len(pos) < min_genes:
        return None
    if not isinstance(null, LCCNull):
        null = LCCNull(layer, n_trials=n_trials, null=null, seed=seed)

    observed = lcc_size(layer.adj, pos)
    rand = null.samples(pos)
    mean = float(rand.mean())
    std = float(rand.std(ddof=1))
    z = (observed - mean) / std if std > 0 else np.nan  # undefined, not 0: see docs/CHANGES.md
    return LCCResult(
        group_id=group_id or "",
        layer_id=layer.id,
        n_genes_in_layer=len(pos),
        lcc_size=observed,
        rand_mean=mean,
        rand_std=std,
        z_score=float(z),
        p_norm=float(norm.sf(z)) if np.isfinite(z) else np.nan,
        p_emp=float(((rand >= observed).sum() + 1) / (len(rand) + 1)),
    )


def bh_adjust(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values (NaNs ignored and kept, as R's p.adjust)."""
    p = np.asarray(p, dtype=float)
    q = np.full_like(p, np.nan)
    ok = np.isfinite(p)
    m = ok.sum()
    if m == 0:
        return q
    pv = p[ok]
    order = np.argsort(pv)[::-1]
    ranked = pv[order] * m / np.arange(m, 0, -1)
    adj = np.minimum(1.0, np.minimum.accumulate(ranked))
    out = np.empty(m)
    out[order] = adj
    q[ok] = out
    return q


def modularity_table(
    multiplex: Multiplex,
    groups: Iterable[GeneGroup],
    n_trials: int = 1000,
    min_genes: int = 10,
    min_lcc: int = 5,
    null: str = "uniform",
    p_value: str = "norm",
    alpha: float = 0.05,
    seed: int = 0,
    nulls: dict[str, LCCNull] | None = None,
) -> pd.DataFrame:
    """LCC modularity for every group x layer, with the paper's significance call.

    Args:
        n_trials: random sets per null (paper: 1000).
        min_genes: minimum group genes in a layer to assess it (paper: 10).
        min_lcc: minimum observed LCC size to call a layer significant (paper: 5).
        null: "uniform" (paper) or "degree".
        p_value: "norm" (paper, normal approximation of z) or "empirical".
        alpha: BH threshold (paper: 0.05), applied over all assessed pairs in the table.
        nulls: optional {layer_id: LCCNull} to share null samples across calls.

    Returns:
        DataFrame with one row per assessed (group, layer): LCCResult fields plus
        ``p_value`` (the chosen one), ``q_value`` and ``significant``.
    """
    if p_value not in ("norm", "empirical"):
        raise ValueError("p_value must be 'norm' or 'empirical'")
    groups = list(groups)
    nulls = nulls if nulls is not None else {}
    rows = []
    for layer in multiplex:
        if layer.id not in nulls:
            nulls[layer.id] = LCCNull(layer, n_trials=n_trials, null=null, seed=seed)
        for g in groups:
            res = lcc_modularity(layer, g, min_genes=min_genes, null=nulls[layer.id])
            if res is not None:
                rows.append(asdict(res))
    cols = list(LCCResult.__dataclass_fields__)
    df = pd.DataFrame(rows, columns=cols)
    df["p_value"] = df["p_norm"] if p_value == "norm" else df["p_emp"]
    df["q_value"] = bh_adjust(df["p_value"].to_numpy())
    df["significant"] = (df["q_value"] < alpha) & (df["lcc_size"] >= min_lcc)
    return df


def significant_layers(table: pd.DataFrame, group_id: str) -> dict[str, float]:
    """{layer_id: z} of the layers called significant for a group (paper's selection)."""
    sub = table[(table["group_id"] == group_id) & table["significant"]]
    return dict(zip(sub["layer_id"], sub["z_score"], strict=True))
