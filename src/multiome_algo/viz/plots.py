"""Core plots: per-group/per-layer modularity heatmap and retrieval ROC."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # headless
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from sklearn.metrics import roc_curve  # noqa: E402

from multiome_algo.lcc import LCCResult  # noqa: E402


def plot_modularity_heatmap(
    results: list[LCCResult],
    out: str | Path | None = None,
    value: str = "z_score",
):
    """Heatmap of LCC modularity (groups x layers).

    Args:
        results: LCCResult objects across groups and layers.
        out: if given, save the figure to this path.
        value: which field to display ("z_score" or "lcc_size").
    """
    groups = sorted({r.group_id for r in results})
    layers = sorted({r.layer_id for r in results})
    gi = {g: i for i, g in enumerate(groups)}
    li = {lid: j for j, lid in enumerate(layers)}

    mat = np.full((len(groups), len(layers)), np.nan)
    for r in results:
        mat[gi[r.group_id], li[r.layer_id]] = getattr(r, value)

    fig, ax = plt.subplots(figsize=(max(6, len(layers) * 0.4), max(4, len(groups) * 0.3)))
    im = ax.imshow(mat, aspect="auto", cmap="RdBu_r")
    ax.set_xticks(range(len(layers)))
    ax.set_xticklabels(layers, rotation=90, fontsize=7)
    ax.set_yticks(range(len(groups)))
    ax.set_yticklabels(groups, fontsize=7)
    ax.set_title(f"LCC modularity ({value})")
    fig.colorbar(im, ax=ax, label=value)
    fig.tight_layout()
    if out:
        fig.savefig(out, dpi=150)
    return fig


def plot_roc(
    labels,
    scores,
    out: str | Path | None = None,
    label: str = "",
):
    """ROC curve for retrieval results (concatenated labels/scores across folds)."""
    fpr, tpr, _ = roc_curve(labels, scores)
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.plot(fpr, tpr, label=label or "RWR")
    ax.plot([0, 1], [0, 1], "--", color="grey", linewidth=0.8)
    ax.set_xlabel("False positive rate")
    ax.set_ylabel("True positive rate")
    ax.set_title("Gene retrieval ROC")
    ax.legend(loc="lower right")
    fig.tight_layout()
    if out:
        fig.savefig(out, dpi=150)
    return fig
