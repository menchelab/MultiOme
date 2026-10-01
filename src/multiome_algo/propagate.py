"""Informed multiplex random walk with restart (RWR).

Port of ``RWR_transitional_matrix.R``, ``RWR.R`` and ``weighted_multiplex_propagation.R``
(Buphamalai et al. 2021).

Supra-transition matrix (paper)
    Every gene has one copy per layer (N genes x L layers states). Block (i, j) - moving
    from layer j to layer i - is

        S_ij = P[i, j] * M_i   if i == j   (M_i: column-normalised adjacency of layer i)
        S_ij = P[i, j] * I     otherwise   (jump to the same gene's copy in layer i)

    and S is then column-normalised. P is the layer transition matrix of `pmat_cal`, a
    Metropolis-Hastings rule whose stationary distribution is proportional to the layer
    weights w (the LCC z-scores of the significant layers):

        P[i, j] = min(1, w_i / w_j) / L   (i != j),   P[j, j] = 1 - sum_{i != j} P[i, j]

    A gene absent from layer j has an empty column in M_j; the final column normalisation
    then routes all its mass to its copies in other layers ("pass-through" copies).

Option ``coupling="present"``: only (gene, layer) pairs where the gene exists in the layer
are states; jumps go only to layers containing the gene; restart only seeds present
copies. This removes pass-through states (which bend the intended weighting) and shrinks
the state space.

Walk:  p_{t+1} = (1 - r) S p_t + r p_0, restart r = 0.7, p_0 = seeds (equal weight, or
given seed weights) replicated across layer copies and normalised to sum 1.

Gene score: per-layer visiting probabilities combined by arithmetic mean (paper; for
"present" coupling absent copies count as 0), geometric mean, or geometric mean of ranks.

Implementation note: S is never materialised. With X the N x L matrix of layer-copy
probabilities, S @ vec(X) is computed blockwise from the L sparse layer matrices and the
L x L matrix P - exactly equal to the paper's matrix, at a fraction of the memory.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import rankdata

from multiome_core.schema import GeneGroup, Multiplex

WEIGHTINGS = ("z", "uniform", "softmax")
COUPLINGS = ("paper", "present")
COMBINES = ("mean", "geomean", "rank_geomean")


def pmat_paper(weights: Iterable[float]) -> np.ndarray:
    """Layer transition matrix of the paper (``pmat_cal``). Columns sum to 1.

    Requires strictly positive weights (the paper only feeds z-scores of significant
    layers, so z > 0 always holds there).
    """
    w = np.asarray(list(weights), dtype=float)
    if np.any(~np.isfinite(w)) or np.any(w <= 0):
        raise ValueError(f"Layer weights must be finite and > 0 (got {w}).")
    L = len(w)
    P = np.minimum(1.0, w[:, None] / w[None, :]) / L
    np.fill_diagonal(P, 0.0)
    np.fill_diagonal(P, 1.0 - P.sum(axis=0))
    return P


def layer_weights(zscores: Mapping[str, float], weighting: str = "z") -> dict[str, float]:
    """Turn per-layer z-scores into the weights fed to `pmat_paper`.

    "z" (paper): raw z-scores.  "uniform": all 1 (the paper's unweighted runs).
    "softmax": exp(z - max z) - an option; much more peaked than the paper's rule.
    """
    ids = list(zscores)
    z = np.array([zscores[i] for i in ids], dtype=float)
    if weighting == "z":
        w = z
    elif weighting == "uniform":
        w = np.ones_like(z)
    elif weighting == "softmax":
        w = np.exp(z - z.max())
    else:
        raise ValueError(f"weighting must be one of {WEIGHTINGS}")
    return dict(zip(ids, w.astype(float), strict=True))


def _column_normalize(a: sp.csr_matrix) -> sp.csr_matrix:
    a = a.tocsc(copy=True).astype(float)
    s = np.asarray(a.sum(axis=0)).ravel()
    inv = np.divide(1.0, s, out=np.zeros_like(s), where=s > 0)
    return (a @ sp.diags(inv)).tocsr()


def normalized_adjacency(
    multiplex: Multiplex, layer_id: str, use_edge_weights: bool = True
) -> sp.csr_matrix:
    """Column-normalised layer adjacency on the multiplex universe (cached on multiplex)."""
    cache = multiplex.__dict__.setdefault("_normalized", {})
    key = (layer_id, use_edge_weights and multiplex[layer_id].weighted)
    if key not in cache:
        a = multiplex.aligned_adjacency(layer_id)
        if not key[1]:
            a = a.astype(bool).astype(float)
        cache[key] = _column_normalize(a)
    return cache[key]


class SupraOperator:
    """Matrix-free supra-transition operator over a set of layers.

    Args:
        multiplex: provides the layers and the shared gene universe.
        layer_ids: the layers to walk (order defines columns of the state matrix).
        weights: {layer_id: weight > 0}; uniform when None.
        coupling: "paper" or "present".
        use_edge_weights: use edge weights of weighted layers (paper layers are binary).
        P: explicit L x L layer transition matrix (overrides `weights`).
    """

    def __init__(
        self,
        multiplex: Multiplex,
        layer_ids: list[str],
        weights: Mapping[str, float] | None = None,
        coupling: str = "paper",
        use_edge_weights: bool = True,
        P: np.ndarray | None = None,
    ):
        if coupling not in COUPLINGS:
            raise ValueError(f"coupling must be one of {COUPLINGS}")
        if not layer_ids:
            raise ValueError("Need at least one layer.")
        self.layer_ids = list(layer_ids)
        self.coupling = coupling
        L = len(self.layer_ids)
        if P is None:
            w = [1.0] * L if weights is None else [weights[lid] for lid in self.layer_ids]
            P = pmat_paper(w)
        self.P = np.asarray(P, dtype=float)

        # State space = the full multiplex universe. Genes absent from every walked layer
        # (inactive) can never receive mass: seeds are restricted to active genes, and
        # results are reported for active genes only.
        self.genes = multiplex.universe
        self.present = multiplex.presence(self.layer_ids)
        self.active = self.present.any(axis=1)
        self.M = [
            normalized_adjacency(multiplex, lid, use_edge_weights) for lid in self.layer_ids
        ]

        d = np.diag(self.P)
        pres = self.present.astype(float)
        mask = pres if coupling == "present" else np.ones_like(pres)
        # column sums of the un-normalised supra matrix, per state (gene, layer j):
        #   sum_{i != j} P[i,j] * mask[g,i]  +  P[j,j] * present[g,j]
        colsum = mask @ self.P - mask * d + pres * d
        self.scale = np.divide(1.0, colsum, out=np.zeros_like(colsum), where=colsum > 0)
        if coupling == "present":
            self.scale *= pres
        self.mask = mask
        self.index = pd.Index(self.genes)

    @property
    def shape(self) -> tuple[int, int]:
        return len(self.genes), len(self.layer_ids)

    def apply(self, X: np.ndarray) -> np.ndarray:
        """One step: S @ vec(X), with X of shape (N, L)."""
        Xs = X * self.scale
        Y = (Xs @ self.P.T) * self.mask  # inter-layer jumps (incl. diagonal, removed below)
        d = np.diag(self.P)
        Y -= Xs * d * self.mask
        for i, M in enumerate(self.M):
            Y[:, i] += d[i] * (M @ Xs[:, i])
        return Y

    def restart_vector(
        self, seeds: Iterable[str] | Mapping[str, float]
    ) -> np.ndarray:
        """p_0 of the paper: seed genes replicated across layer copies, sums to 1."""
        if isinstance(seeds, Mapping):
            names, vals = list(seeds), np.array(list(seeds.values()), dtype=float)
        else:
            names = list(seeds)
            vals = np.ones(len(names))
        pos = self.index.get_indexer(names)
        ok = pos >= 0
        ok[ok] = self.active[pos[ok]]
        if not ok.any():
            raise ValueError("No seed genes are present in the walked layers.")
        X0 = np.zeros(self.shape)
        X0[pos[ok], :] = vals[ok, None]
        if self.coupling == "present":
            X0 *= self.present
        return X0 / X0.sum()

    def walk(
        self,
        seeds: Iterable[str] | Mapping[str, float],
        r: float = 0.7,
        tol: float = 1e-10,
        max_iter: int = 1000,
    ) -> np.ndarray:
        """Stationary RWR probabilities (N x L) for the given seeds."""
        X0 = self.restart_vector(seeds)
        X = X0.copy()
        for _ in range(max_iter):
            X_next = (1.0 - r) * self.apply(X) + r * X0
            if np.abs(X_next - X).sum() < tol:
                return X_next
            X = X_next
        return X


@dataclass
class RankResult:
    """Result of one propagation run.

    Attributes:
        layer_probs: genes x layers visiting probabilities.
        table: per gene - score, rank (1 = best, among non-seeds if seeds removed), is_seed.
        seeds: seed genes used (those present in the walked layers).
        layer_weights: weights used per layer.
    """

    layer_probs: pd.DataFrame
    table: pd.DataFrame
    seeds: frozenset[str]
    layer_weights: dict[str, float]

    def top(self, n: int = 20) -> pd.DataFrame:
        return self.table.nsmallest(n, "rank")


def combine_layers(X: np.ndarray, how: str = "mean") -> np.ndarray:
    """Combine per-layer probabilities (N x L) into one score per gene (higher = better)."""
    if how == "mean":
        return X.mean(axis=1)
    if how == "geomean":
        with np.errstate(divide="ignore"):
            return np.exp(np.log(X).mean(axis=1))
    if how == "rank_geomean":  # R "RankAllAvg": rank all copies jointly, geo-mean ranks
        ranks = rankdata(-X, method="average").reshape(X.shape)
        return -np.exp(np.log(ranks).mean(axis=1))
    raise ValueError(f"combine must be one of {COMBINES}")


def informed_rwr(
    multiplex: Multiplex,
    seeds: Iterable[str] | Mapping[str, float] | GeneGroup,
    layer_weights: Mapping[str, float] | None = None,
    layer_ids: list[str] | None = None,
    coupling: str = "paper",
    combine: str = "mean",
    r: float = 0.7,
    remove_seeds: bool = True,
    operator: SupraOperator | None = None,
    tol: float = 1e-10,
    max_iter: int = 1000,
) -> RankResult:
    """Run (informed) multiplex RWR and rank genes.

    Args:
        multiplex: the layers.
        seeds: seed genes, {gene: weight}, or a GeneGroup.
        layer_weights: {layer_id: weight > 0} (e.g. z-scores of significant layers, see
            `lcc.significant_layers`); its keys define the walked layers unless
            `layer_ids` is given. Uniform over all layers when None.
        coupling: "paper" or "present" (see module docstring).
        combine: "mean" (paper), "geomean" or "rank_geomean".
        r: restart probability (paper: 0.7).
        remove_seeds: rank only non-seed genes (paper: True); seeds get rank NaN.
        operator: a prebuilt SupraOperator to reuse across seed sets (fast CV).
    """
    if isinstance(seeds, GeneGroup):
        seeds = seeds.genes
    if operator is None:
        if layer_ids is None:
            layer_ids = list(layer_weights) if layer_weights else multiplex.layer_ids
        operator = SupraOperator(multiplex, layer_ids, layer_weights, coupling=coupling)
    X = operator.walk(seeds, r=r, tol=tol, max_iter=max_iter)[operator.active]
    genes = operator.genes[operator.active]
    score = combine_layers(X, combine)

    is_seed = pd.Index(genes).isin(list(seeds))
    table = pd.DataFrame({"gene": genes, "score": score, "is_seed": is_seed})
    ranked = table[~table["is_seed"]] if remove_seeds else table
    table["rank"] = np.nan
    table.loc[ranked.index, "rank"] = rankdata(-ranked["score"].to_numpy(), method="average")
    table = table.sort_values(["rank", "score"], ascending=[True, False], na_position="last")

    L = len(operator.layer_ids)
    weights = (
        dict(layer_weights) if layer_weights else dict.fromkeys(operator.layer_ids, 1.0 / L)
    )
    return RankResult(
        layer_probs=pd.DataFrame(X, index=genes, columns=operator.layer_ids),
        table=table.reset_index(drop=True),
        seeds=frozenset(genes[is_seed]),
        layer_weights=weights,
    )
