"""Informed multiplex random walk with restart (RWR).

Ports the core method of Buphamalai et al. 2021: a random walk with restart over a
supra-adjacency matrix built from multiple layers, where each layer's contribution is
weighted by its per-gene-group relevance (a softmax over LCC modularity z-scores).

Construction (matching the paper's "informed" propagation):
  * Each gene exists as one node per layer it appears in -> the supra-state space.
  * Intra-layer edges come from each layer's adjacency, column-normalized (random-walk
    transition), then scaled by that layer's weight.
  * Inter-layer coupling links copies of the same gene across layers, so the walker can
    switch layers; coupling strength is controlled by `delta`.
  * Restart returns to the seed genes (spread across their layer copies) with prob `r`.

Result: a per-gene score = visiting probability averaged over that gene's layer copies.
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp

from multiome_core.schema import GeneGroup, Multiplex


def softmax_weights(zscores: dict[str, float]) -> dict[str, float]:
    """Softmax over layer z-scores -> layer weights summing to 1 (the paper's pi_dm)."""
    ids = list(zscores)
    z = np.array([zscores[i] for i in ids], dtype=float)
    z = z - z.max()  # numerical stability
    e = np.exp(z)
    w = e / e.sum()
    return dict(zip(ids, w.astype(float), strict=True))


def _column_normalize(a: sp.csr_matrix) -> sp.csr_matrix:
    """Column-stochastic normalization (each column sums to 1, empty cols left at 0)."""
    a = a.tocsc(copy=True).astype(float)
    col_sums = np.asarray(a.sum(axis=0)).ravel()
    inv = np.divide(1.0, col_sums, out=np.zeros_like(col_sums), where=col_sums > 0)
    return (a @ sp.diags(inv)).tocsr()


def informed_rwr(
    multiplex: Multiplex,
    seeds: set[str] | GeneGroup,
    layer_weights: dict[str, float] | None = None,
    r: float = 0.7,
    delta: float = 0.5,
    tol: float = 1e-10,
    max_iter: int = 1000,
) -> dict[str, float]:
    """Run informed multiplex RWR and return a {gene: score} ranking over the universe.

    Args:
        multiplex: the layers to walk.
        seeds: seed genes (or a GeneGroup).
        layer_weights: per-layer weights (e.g. softmax of LCC z-scores). Defaults to uniform.
        r: restart probability.
        delta: inter-layer coupling weight (0 = layers never switch, 1 = strong coupling).
        tol: L1 convergence tolerance.
        max_iter: iteration cap.
    """
    if isinstance(seeds, GeneGroup):
        seeds = seeds.genes

    layer_ids = multiplex.layer_ids
    if not layer_ids:
        raise ValueError("Multiplex has no layers.")
    if layer_weights is None:
        layer_weights = {lid: 1.0 / len(layer_ids) for lid in layer_ids}

    # Shared node universe and index.
    universe = sorted(multiplex.universe())
    idx = {g: i for i, g in enumerate(universe)}
    n = len(universe)
    L = len(layer_ids)

    # Per-layer column-normalized transition matrices over the shared universe.
    intra_blocks: list[sp.csr_matrix] = []
    for lid in layer_ids:
        layer = multiplex[lid]
        rows, cols = [], []
        for u, v in layer.graph.edges():
            iu, iv = idx[u], idx[v]
            rows += [iu, iv]
            cols += [iv, iu]  # undirected
        a = sp.csr_matrix(
            (np.ones(len(rows)), (rows, cols)), shape=(n, n)
        )
        intra_blocks.append(_column_normalize(a))

    # Supra transition matrix S of shape (L*n, L*n).
    # Intra-layer (within a layer block): (1-delta) * weight_l * T_l
    # Inter-layer (gene copy across layers): delta coupling, uniform across the other layers.
    w = np.array([layer_weights.get(lid, 0.0) for lid in layer_ids], dtype=float)

    blocks: list[list[sp.spmatrix | None]] = [[None] * L for _ in range(L)]
    inter_scale = delta / (L - 1) if L > 1 else 0.0
    eye = sp.identity(n, format="csr")
    for a in range(L):
        for b in range(L):
            if a == b:
                blocks[a][b] = (1.0 - delta) * w[a] * intra_blocks[a]
            else:
                # walker in layer b can jump to the same gene's copy in layer a
                blocks[a][b] = inter_scale * eye
    S = sp.bmat(blocks, format="csr")
    S = _column_normalize(S)  # ensure column-stochastic over the supra space

    # Restart vector: seeds present in the universe, spread uniformly across layer copies.
    seed_genes = [g for g in seeds if g in idx]
    if not seed_genes:
        raise ValueError("No seed genes are present in the multiplex universe.")
    p0 = np.zeros(L * n)
    for g in seed_genes:
        for a in range(L):
            p0[a * n + idx[g]] = 1.0
    p0 /= p0.sum()

    # Power iteration.
    p = p0.copy()
    for _ in range(max_iter):
        p_next = (1.0 - r) * (S @ p) + r * p0
        if np.abs(p_next - p).sum() < tol:
            p = p_next
            break
        p = p_next

    # Aggregate per gene: average visiting probability across its layer copies.
    p_mat = p.reshape(L, n)
    gene_scores = p_mat.mean(axis=0)
    return {universe[i]: float(gene_scores[i]) for i in range(n)}
