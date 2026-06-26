"""CLI for the algorithm unit.

Subcommands:
  modularity   LCC modularity for gene groups across layers -> TSV (+ optional heatmap)
  cv           k-fold retrieval cross-validation for a gene group -> AUROC/top-k
  rank         informed multiplex RWR for one gene group -> ranked gene scores TSV
"""

from __future__ import annotations

import argparse
import csv
import sys

import numpy as np

from multiome_algo.crossval import retrieval_cv
from multiome_algo.groups import load_orphanet_groups
from multiome_algo.lcc import lcc_modularity
from multiome_algo.propagate import informed_rwr, softmax_weights
from multiome_core.bundle import read_bundle
from multiome_core.legacy import read_multiplex


def _load_multiplex(args):
    if args.bundle:
        return read_bundle(args.bundle, layer_ids=args.layers)
    return read_multiplex(args.edgelists, layer_ids=args.layers)


def _add_source_args(p):
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--bundle", help="Path to a .mpx bundle directory")
    src.add_argument("--edgelists", help="Directory of legacy TSV edge lists")
    p.add_argument("--layers", nargs="*", help="Restrict to these layer ids")
    p.add_argument("--groups", required=True, help="OrphaNet-style gene-group TSV")


def _cmd_modularity(args):
    mpx = _load_multiplex(args)
    groups = load_orphanet_groups(args.groups)
    rng = np.random.default_rng(args.seed)
    results = []
    for g in groups:
        for lid in mpx.layer_ids:
            res = lcc_modularity(mpx[lid], g, n_trials=args.trials, rng=rng)
            if res:
                results.append(res)

    w = csv.writer(sys.stdout, delimiter="\t")
    w.writerow(["group_id", "layer_id", "n_genes_in_layer", "lcc_size", "z_score", "p_value"])
    for r in results:
        w.writerow(
            [r.group_id, r.layer_id, r.n_genes_in_layer, r.lcc_size,
             f"{r.z_score:.4f}", f"{r.p_value:.4g}"]
        )
    if args.heatmap:
        from multiome_algo.viz import plot_modularity_heatmap

        plot_modularity_heatmap(results, out=args.heatmap)


def _cmd_cv(args):
    mpx = _load_multiplex(args)
    groups = {g.id: g for g in load_orphanet_groups(args.groups)}
    group = groups[args.group] if args.group else max(groups.values(), key=len)
    rng = np.random.default_rng(args.seed)

    w = csv.writer(sys.stdout, delimiter="\t")
    w.writerow(["group_id", "weighting", "auroc", "top5", "top10", "top20"])
    for weighting in args.weighting:
        res = retrieval_cv(
            mpx, group, k_folds=args.folds, weighting=weighting,
            r=args.restart, lcc_trials=args.trials, rng=rng,
        )
        w.writerow(
            [res.group_id, res.weighting, f"{res.auroc:.4f}",
             f"{res.top_k.get(5, float('nan')):.3f}",
             f"{res.top_k.get(10, float('nan')):.3f}",
             f"{res.top_k.get(20, float('nan')):.3f}"]
        )


def _cmd_rank(args):
    mpx = _load_multiplex(args)
    groups = {g.id: g for g in load_orphanet_groups(args.groups)}
    group = groups[args.group] if args.group else max(groups.values(), key=len)
    rng = np.random.default_rng(args.seed)

    zs = {}
    for lid in mpx.layer_ids:
        res = lcc_modularity(mpx[lid], group, n_trials=args.trials, rng=rng)
        zs[lid] = res.z_score if res else 0.0
    weights = softmax_weights(zs)
    scores = informed_rwr(mpx, group.genes, layer_weights=weights, r=args.restart)

    w = csv.writer(sys.stdout, delimiter="\t")
    w.writerow(["gene", "score", "is_seed"])
    for gene in sorted(scores, key=scores.get, reverse=True):
        w.writerow([gene, f"{scores[gene]:.6g}", int(gene in group.genes)])


def main(argv=None) -> None:
    p = argparse.ArgumentParser(prog="multiome-algo", description=__doc__)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--trials", type=int, default=200, help="LCC randomization trials")
    sub = p.add_subparsers(dest="cmd", required=True)

    pm = sub.add_parser("modularity", help="LCC modularity across layers")
    _add_source_args(pm)
    pm.add_argument("--heatmap", help="Save modularity heatmap to this path")
    pm.set_defaults(func=_cmd_modularity)

    pc = sub.add_parser("cv", help="k-fold retrieval cross-validation")
    _add_source_args(pc)
    pc.add_argument("--group", help="Group id (default: largest)")
    pc.add_argument("--folds", type=int, default=5)
    pc.add_argument("--restart", type=float, default=0.7)
    pc.add_argument("--weighting", nargs="+", default=["informed", "uniform"])
    pc.set_defaults(func=_cmd_cv)

    pr = sub.add_parser("rank", help="Informed multiplex RWR ranking")
    _add_source_args(pr)
    pr.add_argument("--group", help="Group id (default: largest)")
    pr.add_argument("--restart", type=float, default=0.7)
    pr.set_defaults(func=_cmd_rank)

    args = p.parse_args(argv)
    args.func(args)
