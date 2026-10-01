"""Static figures (matplotlib) for the four core outputs.

    plot_modularity_heatmap   groups x layers LCC z-scores, significant cells marked
    plot_layer_weights        one group's layer relevance and the layers the walk uses
    plot_cv_performance       ROC + per-group AUROC across configurations (paired)
    plot_candidate_ranking    top-ranked candidates with per-layer contributions

Every function returns a Figure and saves it when `out=` is given. Colours: categorical
slots in fixed order (validated for colour-vision deficiency, all pairs, first three
slots), diverging blue-grey-red for signed z-scores, single-hue blue for magnitudes.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # headless
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7",
          "#e34948"]
MARKERS = ["o", "s", "D", "^", "v", "P", "X", "*"]  # secondary encoding beyond 3 series
INK, INK2, MUTED, GRID, AXIS = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
NEUTRAL = "#c3c2b7"
DIVERGING = LinearSegmentedColormap.from_list(
    "multiome_div", ["#104281", "#3987e5", "#f0efec", "#e66767", "#a32b2b"]
)
SEQUENTIAL = LinearSegmentedColormap.from_list(
    "multiome_seq", ["#f0efec", "#cde2fb", "#6da7ec", "#256abf", "#0d366b"]
)


def _style(ax, grid_axis: str = "x") -> None:
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(AXIS)
    ax.tick_params(colors=INK2, labelsize=8, length=0)
    if grid_axis:
        ax.grid(axis=grid_axis, color=GRID, linewidth=0.6)
        ax.set_axisbelow(True)
    ax.xaxis.label.set_color(INK2)
    ax.yaxis.label.set_color(INK2)


def _finish(fig, out):
    fig.tight_layout()
    if out:
        Path(out).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, dpi=200, bbox_inches="tight", facecolor="white")
    return fig


def _series_colors(names: Sequence[str]) -> dict[str, str]:
    if len(names) > len(SERIES):
        raise ValueError(f"At most {len(SERIES)} configurations per plot; facet instead.")
    return dict(zip(names, SERIES, strict=False))


# --------------------------------------------------------------------------- modularity


def plot_modularity_heatmap(
    table: pd.DataFrame,
    layer_tags: Mapping[str, str] | None = None,
    value: str = "z_score",
    clip: float | None = None,
    out: str | Path | None = None,
    labels: Mapping[str, str] | None = None,
):
    """Heatmap of LCC modularity (rows = groups, columns = layers).

    Args:
        table: `modularity_table` output.
        layer_tags: {layer_id: tag} to group columns (e.g. from `Multiplex.summary()`).
        value: column to show ("z_score" by default).
        clip: symmetric colour limit; defaults to the 95th percentile of |value|.
        labels: {group_id: display label}.
    Significant cells (``significant`` column) are marked with a dot; grey = not assessed.
    """
    mat = table.pivot(index="group_id", columns="layer_id", values=value)
    sig = table.pivot(index="group_id", columns="layer_id", values="significant")
    sig = sig.reindex_like(mat).fillna(False).astype(bool)

    tags = layer_tags or {}
    cols = sorted(mat.columns, key=lambda c: (tags.get(c, "~"), -np.nanmean(mat[c]), c))
    rows = mat.loc[:, cols].mean(axis=1).sort_values(ascending=False).index
    mat, sig = mat.loc[rows, cols], sig.loc[rows, cols]

    finite = np.abs(mat.to_numpy()[np.isfinite(mat.to_numpy())])
    lim = clip or (float(np.percentile(finite, 95)) if finite.size else 1.0)
    lim = max(lim, 1e-6)
    fig, ax = plt.subplots(figsize=(max(6, 0.22 * len(cols) + 3), max(3, 0.24 * len(rows) + 1.5)))
    cmap = DIVERGING.with_extremes(bad="#f5f5f3")
    im = ax.imshow(np.ma.masked_invalid(mat.to_numpy()), aspect="auto", cmap=cmap,
                   norm=TwoSlopeNorm(vcenter=0, vmin=-lim, vmax=lim))
    yy, xx = np.nonzero(sig.to_numpy())
    ax.scatter(xx, yy, s=6, color=INK, linewidths=0, label="significant (BH q < 0.05)")

    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(cols, rotation=90, fontsize=7)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([(labels or {}).get(r, r) for r in rows], fontsize=7)
    _style(ax, grid_axis="")
    # separators between tag groups
    if tags:
        group_of = [tags.get(c, "") for c in cols]
        for k in range(1, len(cols)):
            if group_of[k] != group_of[k - 1]:
                ax.axvline(k - 0.5, color="white", linewidth=2)
    cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01, extend="both")
    cb.set_label("LCC z-score" if value == "z_score" else value, color=INK2, fontsize=8)
    cb.outline.set_visible(False)
    cb.ax.tick_params(labelsize=7, colors=INK2, length=0)
    ax.legend(loc="lower right", bbox_to_anchor=(1, 1), frameon=False, fontsize=7,
              labelcolor=INK2, handletextpad=0.2, borderaxespad=0.2)
    ax.set_title("Gene-group modularity per layer", loc="left", color=INK, fontsize=10)
    return _finish(fig, out)


def plot_layer_weights(
    table: pd.DataFrame,
    group_id: str,
    out: str | Path | None = None,
):
    """Lollipop of one group's LCC z-score per layer. Significant layers (those the
    informed walk uses) are coloured and annotated with their share of the stationary
    layer distribution (proportional to z under the paper's rule)."""
    sub = table[table["group_id"] == group_id].sort_values("z_score")
    sub = sub[np.isfinite(sub["z_score"])]
    sel = sub["significant"].to_numpy()
    share = np.where(sel, sub["z_score"], 0.0)
    share = share / share.sum() if share.sum() > 0 else share

    fig, ax = plt.subplots(figsize=(6, max(2.5, 0.2 * len(sub) + 1)))
    y = np.arange(len(sub))
    colors = np.where(sel, SERIES[0], NEUTRAL)
    ax.hlines(y, 0, sub["z_score"], color=colors, linewidth=2)
    ax.scatter(sub["z_score"], y, color=colors, s=36, zorder=3, edgecolor="white", linewidth=1.5)
    ax.axvline(0, color=AXIS, linewidth=0.8)
    for yi, z, s, on in zip(y, sub["z_score"], share, sel, strict=True):
        if on:
            ax.text(z, yi, f"  {s:.0%}", va="center", fontsize=7, color=INK2)
    ax.set_yticks(y)
    ax.set_yticklabels(sub["layer_id"], fontsize=7)
    ax.set_xlabel("LCC z-score")
    _style(ax, "x")
    ax.scatter([], [], color=SERIES[0], s=36, label="used by informed walk (share of walk)")
    ax.scatter([], [], color=NEUTRAL, s=36, label="not significant")
    ax.legend(frameon=False, fontsize=7, loc="lower right", labelcolor=INK2)
    ax.set_title(f"Layer relevance - {group_id}", loc="left", color=INK, fontsize=10)
    return _finish(fig, out)


# --------------------------------------------------------------------------- CV


def plot_cv_performance(cv, out: str | Path | None = None, configs: Sequence[str] | None = None):
    """Two panels: (a) ROC averaged over groups per configuration; (b) per-group median
    AUROC across folds, one line per group connecting the configurations (paired)."""
    summary = cv.summary()
    configs = list(configs or dict.fromkeys(summary["config"]))
    colors = _series_colors(configs)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.2),
                                   gridspec_kw={"width_ratios": [1, 1.3]})
    for k, c in enumerate(configs):
        roc = cv.roc[cv.roc["config"] == c].groupby("fpr")["tpr"].mean()
        auc = summary.loc[summary["config"] == c, "auroc_median"].median()
        ax1.plot(roc.index, roc.values, color=colors[c], linewidth=2,
                 label=f"{c} (median AUROC {auc:.2f})",
                 marker=MARKERS[k] if len(configs) > 3 else None, markevery=10, markersize=5)
    ax1.plot([0, 1], [0, 1], color=AXIS, linewidth=0.8, linestyle="--")
    ax1.set_xlabel("False positive rate")
    ax1.set_ylabel("True positive rate")
    ax1.set_xlim(0, 1)
    ax1.set_ylim(0, 1.01)
    _style(ax1, "both")
    ax1.legend(frameon=False, fontsize=7, loc="lower right", labelcolor=INK2)
    ax1.set_title("a  Retrieval ROC (mean over groups)", loc="left", color=INK, fontsize=10)

    wide = summary.pivot(index="group_id", columns="config", values="auroc_median")
    wide = wide[[c for c in configs if c in wide]]
    x = np.arange(wide.shape[1])
    # one horizontal offset per group, shared by its points and its connecting line
    jitter = (np.random.default_rng(0).random(len(wide)) - 0.5) * 0.12
    xs = x[None, :] + jitter[:, None]
    for gi in range(len(wide)):
        ax2.plot(xs[gi], wide.iloc[gi].to_numpy(), color=GRID, linewidth=1, zorder=1)
    for k, c in enumerate(wide.columns):
        ok = wide[c].notna().to_numpy()
        vals = wide[c].to_numpy()[ok]
        ax2.scatter(xs[ok, k], vals, color=colors[c], s=24, zorder=3,
                    marker=MARKERS[k] if len(configs) > 3 else "o",
                    edgecolor="white", linewidth=1)
        med = np.median(vals) if len(vals) else np.nan
        ax2.hlines(med, k - 0.25, k + 0.25, color=INK, linewidth=2, zorder=4)
        ax2.text(k + 0.28, med, f"{med:.2f}", va="center", fontsize=7, color=INK2)
    ax2.set_xticks(x)
    ax2.set_xticklabels(wide.columns, fontsize=8)
    ax2.set_ylabel("Median AUROC across folds")
    _style(ax2, "y")
    ax2.set_title(f"b  Per group ({cv.protocol} protocol, n = {len(wide)})", loc="left",
                  color=INK, fontsize=10)
    return _finish(fig, out)


# --------------------------------------------------------------------------- ranking


def plot_candidate_ranking(
    result,
    top_n: int = 30,
    highlight: set[str] | None = None,
    out: str | Path | None = None,
    title: str = "Top candidates",
):
    """Top-`top_n` non-seed genes by score, with each layer's share of the gene's score.

    Args:
        result: `RankResult` from `informed_rwr`.
        highlight: genes to mark (e.g. known/held-out disease genes).
    """
    top = result.table[~result.table["is_seed"]].nsmallest(top_n, "rank")
    probs = result.layer_probs.loc[top["gene"]]
    contrib = probs.div(probs.sum(axis=1).replace(0, np.nan), axis=0)
    order = contrib.sum().sort_values(ascending=False).index
    contrib = contrib[order]
    highlight = highlight or set()

    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(4 + 0.25 * len(order) + 3, max(3, 0.24 * top_n + 1.2)), sharey=True,
        gridspec_kw={"width_ratios": [3, max(1.0, 0.25 * len(order))]},
    )
    y = np.arange(len(top))
    is_hl = top["gene"].isin(highlight).to_numpy()
    ax1.barh(y, top["score"], color=np.where(is_hl, SERIES[1], SERIES[0]), height=0.7)
    ax1.set_yticks(y)
    ax1.set_yticklabels(top["gene"], fontsize=7)
    ax1.invert_yaxis()
    ax1.set_xlabel("Score (mean visiting probability)")
    ax1.ticklabel_format(axis="x", style="sci", scilimits=(-2, 2))
    _style(ax1, "x")
    if highlight:
        ax1.legend(handles=[Patch(color=SERIES[0], label="candidate"),
                            Patch(color=SERIES[1], label="highlighted (known)")],
                   frameon=False, fontsize=7, loc="lower right", labelcolor=INK2)
    ax1.set_title(title, loc="left", color=INK, fontsize=10)

    im = ax2.imshow(contrib.to_numpy(), aspect="auto", cmap=SEQUENTIAL, vmin=0,
                    vmax=max(0.01, float(np.nanmax(contrib.to_numpy()))))
    ax2.set_xticks(range(len(order)))
    ax2.set_xticklabels(order, rotation=90, fontsize=7)
    _style(ax2, "")
    cb = fig.colorbar(im, ax=ax2, fraction=0.05, pad=0.02)
    cb.set_label("Layer share of score", color=INK2, fontsize=8)
    cb.outline.set_visible(False)
    cb.ax.tick_params(labelsize=7, colors=INK2, length=0)
    ax2.set_title("Per-layer contribution", loc="left", color=INK, fontsize=10)
    return _finish(fig, out)
