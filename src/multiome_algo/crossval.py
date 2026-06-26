"""Cross-validated gene retrieval — the paper's performance assessment.

k-fold CV per gene group: hold out a fold of the group's genes, seed propagation with
the rest, and measure how well the held-out genes are recovered among the ranking of
all non-seed genes (AUROC + top-k). Compares informed multiplex weighting against
baselines (uniform layer weights, single layer).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from sklearn.metrics import roc_auc_score

from multiome_algo.lcc import lcc_modularity
from multiome_algo.propagate import informed_rwr, softmax_weights
from multiome_core.schema import GeneGroup, Multiplex


@dataclass
class CVResult:
    group_id: str
    weighting: str
    auroc: float
    top_k: dict[int, float] = field(default_factory=dict)  # k -> fraction held-out in top-k
    n_folds: int = 0
    n_genes_evaluated: int = 0


def _layer_weights_for(
    multiplex: Multiplex,
    seeds: set[str],
    group_id: str,
    weighting: str,
    n_trials: int,
    rng: np.random.Generator,
) -> dict[str, float]:
    """Build layer weights under a given weighting scheme using only the seed genes."""
    if weighting == "uniform":
        ids = multiplex.layer_ids
        return {lid: 1.0 / len(ids) for lid in ids}
    if weighting == "informed":
        seed_group = GeneGroup(id=group_id, genes=seeds)
        zs = {}
        for lid in multiplex.layer_ids:
            res = lcc_modularity(multiplex[lid], seed_group, n_trials=n_trials, rng=rng)
            zs[lid] = res.z_score if res else 0.0
        return softmax_weights(zs)
    raise ValueError(f"Unknown weighting: {weighting!r}")


def _evaluate_fold(
    multiplex: Multiplex,
    seeds: set[str],
    held_out: set[str],
    universe: set[str],
    weights: dict[str, float],
    r: float,
    ks: tuple[int, ...],
) -> tuple[list[float], list[int], dict[int, int]]:
    """Return (scores, labels, top_k_hits) over candidate genes (universe minus seeds)."""
    scores = informed_rwr(multiplex, seeds, layer_weights=weights, r=r)
    candidates = [g for g in universe if g not in seeds]
    y_score = np.array([scores.get(g, 0.0) for g in candidates])
    y_true = np.array([1 if g in held_out else 0 for g in candidates])

    ranked = [g for g in sorted(candidates, key=lambda g: scores.get(g, 0.0), reverse=True)]
    hits = {k: sum(1 for g in ranked[:k] if g in held_out) for k in ks}
    return list(y_score), list(y_true), hits


def retrieval_cv(
    multiplex: Multiplex,
    group: GeneGroup,
    k_folds: int = 5,
    weighting: str = "informed",
    r: float = 0.7,
    ks: tuple[int, ...] = (5, 10, 20),
    lcc_trials: int = 200,
    rng: np.random.Generator | None = None,
) -> CVResult:
    """k-fold retrieval CV for one gene group.

    Pooled AUROC across folds (scores/labels concatenated). top-k reported as the mean
    fraction of held-out genes recovered in the top k per fold.
    """
    rng = rng or np.random.default_rng(0)
    universe = multiplex.universe()
    present = sorted(group.in_universe(universe))
    if len(present) < k_folds:
        raise ValueError(
            f"Group {group.id!r} has only {len(present)} genes in universe; need >= {k_folds}."
        )

    perm = rng.permutation(present)
    folds = np.array_split(perm, k_folds)

    all_scores: list[float] = []
    all_labels: list[int] = []
    topk_fracs: dict[int, list[float]] = {k: [] for k in ks}

    for fold in folds:
        held_out = set(fold.tolist())
        seeds = set(present) - held_out
        if not seeds or not held_out:
            continue
        weights = _layer_weights_for(
            multiplex, seeds, group.id, weighting, lcc_trials, rng
        )
        s, y, hits = _evaluate_fold(multiplex, seeds, held_out, universe, weights, r, ks)
        all_scores += s
        all_labels += y
        for k in ks:
            topk_fracs[k].append(hits[k] / len(held_out))

    auroc = (
        float(roc_auc_score(all_labels, all_scores))
        if len(set(all_labels)) > 1
        else float("nan")
    )
    top_k = {k: float(np.mean(v)) if v else float("nan") for k, v in topk_fracs.items()}
    return CVResult(
        group_id=group.id,
        weighting=weighting,
        auroc=auroc,
        top_k=top_k,
        n_folds=k_folds,
        n_genes_evaluated=len(present),
    )
