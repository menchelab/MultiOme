"""Cross-validated gene retrieval (the paper's performance assessment).

For each gene group, its genes are split into k folds (paper: 10). Each fold is held out
in turn: the remaining genes seed the walk and the held-out genes should rank highly among
all candidates (genes of the walked layers, minus seeds). Per fold we record AUROC (the
paper's metric; reported as median/IQR across folds) and top-k recovery.

Several configurations are compared on the *same* folds (paired), e.g.

    informed     significant layers, weighted by z (the paper's method)
    all_uniform  all layers, unweighted
    ppi          PPI layer alone

Protocols - where layer selection and weights come from:
    "train" (default): recomputed from the training genes of each fold. No information
        from held-out genes leaks into the model.
    "paper": computed once from the full gene group, as in the published analysis. The
        held-out genes then influence which layers are used and how they are weighted,
        which inflates performance; use only to reproduce the paper.
"""

from __future__ import annotations

import warnings
import zlib
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, roc_curve

from multiome_algo.lcc import LCCNull, modularity_table, significant_layers
from multiome_algo.propagate import SupraOperator, combine_layers, layer_weights
from multiome_core.schema import GeneGroup, Multiplex, filter_groups

PROTOCOLS = ("train", "paper")
ROC_GRID = np.linspace(0, 1, 101)


@dataclass
class CVConfig:
    """One propagation configuration to evaluate.

    Attributes:
        name: label in results/plots.
        layers: "significant" (per group, from LCC modularity), "all", or a list of ids.
        weighting: "z" (paper), "uniform" or "softmax"; applied to the selected layers.
            Non-"significant" selections have no z-scores and are always uniform unless
            weighting is "z", in which case the group's z-scores are used where positive.
        coupling: "paper" or "present".
        combine: "mean" (paper), "geomean" or "rank_geomean".
    """

    name: str
    layers: str | list[str] = "significant"
    weighting: str = "z"
    coupling: str = "paper"
    combine: str = "mean"


def paper_configs(multiplex: Multiplex) -> list[CVConfig]:
    """The paper's main comparison: informed multiplex vs all layers vs PPI alone."""
    configs = [
        CVConfig("informed", "significant", "z"),
        CVConfig("all_uniform", "all", "uniform"),
    ]
    if "ppi" in multiplex:
        configs.append(CVConfig("ppi", ["ppi"], "uniform"))
    return configs


@dataclass
class CVResult:
    """Per-fold results of `retrieval_cv`.

    Attributes:
        folds: one row per (group, config, fold): auroc, top-k fractions, sizes, layers.
        held_out: rank of every held-out gene in every fold.
        roc: mean ROC curve per (group, config) on a fixed FPR grid.
        modularity: the LCC tables used (one per fold for "train", one for "paper").
    """

    folds: pd.DataFrame
    held_out: pd.DataFrame
    roc: pd.DataFrame
    modularity: dict[int | str, pd.DataFrame] = field(default_factory=dict)
    protocol: str = "train"

    def summary(self) -> pd.DataFrame:
        """Per group x config: median and IQR of AUROC across folds (paper's Fig. 5b),
        plus mean top-k recovery."""
        topk = [c for c in self.folds.columns if c.startswith("top")]
        g = self.folds.groupby(["group_id", "config"], sort=False)
        out = g["auroc"].agg(
            auroc_median="median",
            auroc_q25=lambda s: s.quantile(0.25),
            auroc_q75=lambda s: s.quantile(0.75),
            auroc_mean="mean",
            n_folds="count",
        )
        return out.join(g[topk].mean()).reset_index()


def make_folds(genes: Iterable[str], k: int, seed: int, group_id: str) -> list[list[str]]:
    """Random, near-equal folds; deterministic per (seed, group id)."""
    genes = sorted(genes)
    rng = np.random.default_rng([seed, zlib.crc32(group_id.encode())])
    perm = rng.permutation(len(genes))
    return [[genes[i] for i in part] for part in np.array_split(perm, k)]


def _select(
    config: CVConfig, multiplex: Multiplex, table: pd.DataFrame | None, group_id: str
) -> dict[str, float]:
    """Layer weights {layer: w > 0} for a config; empty if nothing is selected."""
    if config.layers == "significant":
        zs = significant_layers(table, group_id) if table is not None else {}
        return layer_weights(zs, config.weighting) if zs else {}
    ids = multiplex.layer_ids if config.layers == "all" else list(config.layers)
    if config.weighting == "uniform":
        return dict.fromkeys(ids, 1.0)
    if table is None:
        raise ValueError(f"Config {config.name!r} needs z-scores for weighting.")
    sub = table[table["group_id"] == group_id].set_index("layer_id")["z_score"]
    zs = {lid: float(sub.get(lid, np.nan)) for lid in ids}
    zs = {lid: z for lid, z in zs.items() if np.isfinite(z) and z > 0}
    return layer_weights(zs, config.weighting) if zs else {}


def retrieval_cv(
    multiplex: Multiplex,
    groups: Sequence[GeneGroup],
    configs: Sequence[CVConfig] | None = None,
    k_folds: int = 10,
    protocol: str = "train",
    r: float = 0.7,
    top_k: Sequence[int] = (10, 50, 100),
    min_size: int | None = 20,
    max_size: int | None = 2000,
    seed: int = 0,
    modularity: pd.DataFrame | None = None,
    lcc_kwargs: dict | None = None,
    progress: bool = False,
) -> CVResult:
    """k-fold retrieval cross-validation over gene groups and configurations.

    Args:
        multiplex: all candidate layers.
        groups: gene groups; filtered to [min_size, max_size] total genes (paper: 20-2000).
        configs: configurations to compare on identical folds (default: `paper_configs`).
        k_folds: folds per group (paper: 10).
        protocol: "train" (default, no leakage) or "paper" (selection on the full group).
        r: restart probability.
        top_k: report the fraction of held-out genes ranked within the top k.
        seed: controls folds and LCC nulls.
        modularity: precomputed full-group `modularity_table` (protocol "paper").
        lcc_kwargs: passed to `modularity_table` (n_trials, null, p_value, alpha, ...).
    """
    if protocol not in PROTOCOLS:
        raise ValueError(f"protocol must be one of {PROTOCOLS}")
    configs = list(configs or paper_configs(multiplex))
    lcc_kwargs = {"seed": seed, **(lcc_kwargs or {})}
    universe = multiplex.gene_set()

    groups = filter_groups(groups, min_size, max_size)
    folds_by_group: dict[str, list[list[str]]] = {}
    for g in groups:
        present = g.in_universe(universe)
        if len(present) < k_folds:
            warnings.warn(f"Skipping {g.id!r}: {len(present)} genes < {k_folds} folds.",
                          stacklevel=2)
            continue
        folds_by_group[g.id] = make_folds(present, k_folds, seed, g.id)
    groups = [g for g in groups if g.id in folds_by_group]

    # LCC tables: once on full groups ("paper"), or per fold on training genes ("train");
    # nulls are shared across folds so each (layer, size) is sampled only once.
    nulls: dict[str, LCCNull] = {}
    tables: dict[int | str, pd.DataFrame] = {}
    if protocol == "paper":
        tables["full"] = (
            modularity if modularity is not None
            else modularity_table(multiplex, groups, nulls=nulls, **lcc_kwargs)
        )

    fold_rows, ho_rows, roc_rows = [], [], []
    operators: dict[tuple, SupraOperator] = {}
    for f in range(k_folds):
        if protocol == "train":
            train_groups = [
                GeneGroup(id=g.id, genes=set().union(*(fl for j, fl in
                          enumerate(folds_by_group[g.id]) if j != f)))
                for g in groups
            ]
            tables[f] = modularity_table(multiplex, train_groups, nulls=nulls, **lcc_kwargs)
        table = tables["full"] if protocol == "paper" else tables[f]

        for g in groups:
            held = set(folds_by_group[g.id][f])
            train = set().union(*(fl for j, fl in enumerate(folds_by_group[g.id]) if j != f))
            for cfg in configs:
                weights = _select(cfg, multiplex, table, g.id)
                row = {"group_id": g.id, "config": cfg.name, "fold": f,
                       "n_layers": len(weights), "layers": ";".join(weights)}
                if not weights:
                    fold_rows.append({**row, "auroc": np.nan, "status": "no_layers"})
                    continue
                key = (tuple(weights.items()), cfg.coupling)
                if key not in operators:
                    operators[key] = SupraOperator(
                        multiplex, list(weights), weights, coupling=cfg.coupling
                    )
                op = operators[key]
                active = op.active
                seeds = [s for s in train if s in op.index]
                X = op.walk(seeds, r=r)
                score = combine_layers(X, cfg.combine)

                genes = op.genes
                is_seed = op.index.isin(seeds)
                cand = active & ~is_seed
                labels = op.index.isin(held) & cand
                y, s = labels[cand], score[cand]
                if y.sum() == 0 or y.all():
                    fold_rows.append({**row, "auroc": np.nan, "status": "no_positives"})
                    continue
                order = np.argsort(-s, kind="stable")
                ranks = pd.Series(s).rank(ascending=False, method="average").to_numpy()
                n_pos = int(y.sum())
                row.update(
                    auroc=float(roc_auc_score(y, s)),
                    status="ok",
                    n_seeds=len(seeds),
                    n_held_out=n_pos,
                    n_held_out_missing=len(held) - n_pos,
                    n_candidates=int(cand.sum()),
                    **{f"top{k}": float(y[order[:k]].sum() / n_pos) for k in top_k},
                )
                fold_rows.append(row)
                cand_genes = genes[cand]
                for gi in np.flatnonzero(y):
                    ho_rows.append({"group_id": g.id, "config": cfg.name, "fold": f,
                                    "gene": cand_genes[gi], "rank": ranks[gi],
                                    "n_candidates": row["n_candidates"]})
                fpr, tpr, _ = roc_curve(y, s)
                roc_rows.append((g.id, cfg.name, f, np.interp(ROC_GRID, fpr, tpr)))
        # operators depend on fold-specific weights under "train"; drop to bound memory
        if protocol == "train":
            operators.clear()
        if progress:
            print(f"fold {f + 1}/{k_folds} done", flush=True)

    roc = pd.DataFrame(
        [
            {"group_id": gid, "config": c, "fpr": x, "tpr": t}
            for (gid, c), grp in pd.DataFrame(
                roc_rows, columns=["group_id", "config", "fold", "tpr"]
            ).groupby(["group_id", "config"], sort=False)
            for x, t in zip(ROC_GRID, np.mean(np.vstack(grp["tpr"].to_list()), axis=0),
                            strict=True)
        ]
    ) if roc_rows else pd.DataFrame(columns=["group_id", "config", "fpr", "tpr"])

    return CVResult(
        folds=pd.DataFrame(fold_rows),
        held_out=pd.DataFrame(ho_rows),
        roc=roc,
        modularity=tables,
        protocol=protocol,
    )
