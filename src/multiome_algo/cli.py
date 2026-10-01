"""Command-line interface.

    multiome modularity  LCC modularity of gene groups in every layer (+ heatmap)
    multiome rank        informed multiplex RWR from a gene group or seed list
    multiome cv          k-fold retrieval cross-validation, paired across configurations
    multiome run         any of the above from a YAML file of the same options

Inputs are either your own files (``--network DIR --groups FILE``) or the data shipped
with the 2021 paper (``--paper``). Every command writes TSVs and figures into ``--out``.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import pandas as pd

from multiome_algo.crossval import CVConfig, paper_configs, retrieval_cv
from multiome_algo.lcc import modularity_table, significant_layers
from multiome_algo.propagate import COMBINES, COUPLINGS, WEIGHTINGS, informed_rwr, layer_weights
from multiome_core.io import coverage_report, read_gene_groups, read_multiplex
from multiome_core.schema import GeneGroup, filter_groups


def _log(msg: str) -> None:
    print(f"[multiome] {msg}", file=sys.stderr, flush=True)


# --------------------------------------------------------------------------- inputs


def _add_input_args(p: argparse.ArgumentParser, groups_required: bool = True) -> None:
    src = p.add_argument_group("input")
    s = src.add_mutually_exclusive_group(required=True)
    s.add_argument("--network", help="folder of edge-list files, one per layer")
    s.add_argument("--paper", action="store_true", help="use the data shipped with the paper")
    src.add_argument("--groups", help="gene sets: GMT, long (group, gene) or wide table")
    src.add_argument("--layers", nargs="+", help="only use these layer ids")
    src.add_argument("--edge-weights", action="store_true",
                     help="read a third column as edge weight")
    src.add_argument("--normalize-ids", action=argparse.BooleanOptionalAction, default=None,
                     help="map gene symbols to current HGNC (default: on for --paper only)")
    src.add_argument("--min-size", type=int, default=20, help="minimum group size (paper: 20)")
    src.add_argument("--max-size", type=int, default=2000,
                     help="maximum group size (paper: 2000)")
    src.add_argument("--no-cache", action="store_true", help="do not cache parsed layers")
    p.add_argument("--out", default="multiome_out", help="output folder")
    p.add_argument("--seed", type=int, default=0, help="random seed (folds, nulls)")
    p.set_defaults(_groups_required=groups_required)


def _add_lcc_args(p: argparse.ArgumentParser) -> None:
    g = p.add_argument_group("layer significance (LCC modularity)")
    g.add_argument("--trials", type=int, default=1000, help="random sets per null (paper: 1000)")
    g.add_argument("--null", choices=["uniform", "degree"], default="uniform",
                   help="random-set null (paper: uniform)")
    g.add_argument("--p-value", choices=["norm", "empirical"], default="norm",
                   help="p from normal approx. of z (paper) or empirical")
    g.add_argument("--alpha", type=float, default=0.05, help="BH q threshold (paper: 0.05)")


def _add_walk_args(p: argparse.ArgumentParser) -> None:
    g = p.add_argument_group("propagation")
    g.add_argument("--weighting", choices=WEIGHTINGS, default="z",
                   help="layer weights from z-scores (paper: z)")
    g.add_argument("--coupling", choices=COUPLINGS, default="paper",
                   help="inter-layer coupling (paper: every gene in every layer)")
    g.add_argument("--combine", choices=COMBINES, default="mean",
                   help="combine per-layer probabilities (paper: mean)")
    g.add_argument("--restart", type=float, default=0.7, help="restart probability (paper: 0.7)")


def _load(args):
    """Returns (multiplex, groups) after optional normalisation; writes coverage.tsv."""
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    cache = {} if not args.no_cache else {"cache_dir": None}
    t0 = time.time()
    if args.paper:
        from multiome_core.legacy import read_paper_groups, read_paper_multiplex

        mpx = read_paper_multiplex(layer_ids=args.layers, **cache)
        groups = read_gene_groups(args.groups) if args.groups else read_paper_groups()
        normalize = True if args.normalize_ids is None else args.normalize_ids
    else:
        mpx = read_multiplex(args.network, layer_ids=args.layers, weights=args.edge_weights,
                             **cache)
        groups = read_gene_groups(args.groups) if args.groups else []
        normalize = bool(args.normalize_ids)
    _log(f"loaded {len(mpx)} layers, {len(mpx.universe):,} genes, {len(groups)} groups "
         f"({time.time() - t0:.1f}s)")
    if args._groups_required and not groups:
        raise SystemExit("--groups is required (or use --paper)")

    if normalize:
        from multiome_core.ids import normalize_symbols

        mpx, groups, report = normalize_symbols(mpx, groups)
        report.to_csv(out / "symbol_mapping.tsv", sep="\t", index=False)
        changed = (report["output"] != "").sum()
        _log(f"normalised symbols: {changed} renamed, {len(report) - changed} unresolved "
             f"(see symbol_mapping.tsv)")
    if groups:
        cov = coverage_report(mpx, groups)
        cov.to_csv(out / "coverage.tsv", sep="\t", index=False)
        n_missing = int(cov["n_missing"].sum())
        if n_missing:
            _log(f"{n_missing} group-gene memberships not in any layer (see coverage.tsv)")
    return mpx, groups


def _select_groups(groups: list[GeneGroup], args) -> list[GeneGroup]:
    groups = filter_groups(groups, args.min_size, args.max_size)
    wanted = getattr(args, "group", None)
    if wanted:
        by_key = {g.id: g for g in groups} | {g.label: g for g in groups}
        missing = [w for w in wanted if w not in by_key]
        if missing:
            raise SystemExit(f"unknown group(s) (or outside size filter): {missing}")
        groups = [by_key[w] for w in wanted]
    if not groups:
        raise SystemExit("no gene groups left after the size filter")
    return groups


def _lcc_kwargs(args) -> dict:
    return {"n_trials": args.trials, "null": args.null, "p_value": args.p_value,
            "alpha": args.alpha, "seed": args.seed}


def _tags(mpx) -> dict[str, str]:
    return {lyr.id: (lyr.tags[0] if lyr.tags else "") for lyr in mpx}


# --------------------------------------------------------------------------- commands


def cmd_modularity(args) -> None:
    mpx, groups = _load(args)
    groups = _select_groups(groups, args)
    t0 = time.time()
    table = modularity_table(mpx, groups, **_lcc_kwargs(args))
    out = Path(args.out)
    table.to_csv(out / "modularity.tsv", sep="\t", index=False)
    _log(f"modularity: {len(groups)} groups x {len(mpx)} layers, "
         f"{int(table['significant'].sum())} significant pairs ({time.time() - t0:.1f}s)")
    if not args.no_plots:
        from multiome_algo.viz import plot_layer_weights, plot_modularity_heatmap

        labels = {g.id: g.label for g in groups}
        plot_modularity_heatmap(table, layer_tags=_tags(mpx), labels=labels,
                                out=out / "modularity_heatmap.png")
        import matplotlib.pyplot as plt

        for g in groups if args.per_group_plots else []:
            plt.close(plot_layer_weights(table, g.id, out=out / "layers" / f"{_safe(g.id)}.png"))
    _log(f"wrote {out}/modularity.tsv")


def _seed_group(args, groups) -> GeneGroup:
    if args.seeds:
        p = Path(args.seeds[0])
        genes = p.read_text().split() if len(args.seeds) == 1 and p.exists() else args.seeds
        return GeneGroup(id=args.name or "seeds", genes=genes, source="cli")
    if not args.group or len(args.group) != 1:
        raise SystemExit("rank needs exactly one --group, or --seeds")
    return _select_groups(groups, args)[0]


def cmd_rank(args) -> None:
    mpx, groups = _load(args)
    group = _seed_group(args, groups)
    seeds = group.in_universe(mpx.gene_set())
    _log(f"{group.id}: {len(seeds)}/{len(group)} seed genes in the multiplex")

    if args.layer_selection == "all":
        weights = dict.fromkeys(mpx.layer_ids, 1.0)
        table = None
    else:
        if args.modularity:
            table = pd.read_csv(args.modularity, sep="\t", dtype={"group_id": str})
            if group.id not in set(table["group_id"]):
                raise SystemExit(f"{group.id!r} not in {args.modularity}")
        else:
            # significance is BH-corrected over the table: a single group here is its own
            # family (the paper corrected over all groups at once - pass --modularity)
            table = modularity_table(mpx, [group], **_lcc_kwargs(args))
        zs = significant_layers(table, group.id)
        if not zs:
            raise SystemExit(f"no significant layer for {group.id!r}; try --layer-selection all")
        weights = layer_weights(zs, args.weighting)
    _log("layers: " + ", ".join(f"{k} ({v:.2f})" for k, v in weights.items()))

    res = informed_rwr(mpx, seeds, layer_weights=weights, coupling=args.coupling,
                       combine=args.combine, r=args.restart)
    out = Path(args.out)
    tab = res.table.join(res.layer_probs.add_prefix("p:"), on="gene")
    tab.to_csv(out / f"ranking_{_safe(group.id)}.tsv", sep="\t", index=False)
    if table is not None:
        table.to_csv(out / f"modularity_{_safe(group.id)}.tsv", sep="\t", index=False)
    (out / f"weights_{_safe(group.id)}.json").write_text(json.dumps(weights, indent=2))
    if not args.no_plots:
        from multiome_algo.viz import plot_candidate_ranking, plot_layer_weights

        known = set(Path(args.highlight).read_text().split()) if args.highlight else None
        plot_candidate_ranking(res, top_n=args.top, highlight=known, title=group.label,
                               out=out / f"ranking_{_safe(group.id)}.png")
        if table is not None:
            plot_layer_weights(table, group.id, out=out / f"layers_{_safe(group.id)}.png")
    print(res.top(args.top)[["gene", "score", "rank"]].to_string(index=False))


def _configs(args, mpx) -> list[CVConfig]:
    if not args.configs:
        cfgs = paper_configs(mpx)
        for c in cfgs:
            c.coupling, c.combine = args.coupling, args.combine
            if c.name == "informed":
                c.weighting = args.weighting
        return cfgs
    cfgs = []
    for spec in args.configs:
        # name=layers[:weighting]  layers: significant | all | id+id+...
        name, _, rest = spec.partition("=")
        layers, _, weighting = (rest or name).partition(":")
        sel = layers if layers in ("significant", "all") else layers.split("+")
        default_w = args.weighting if sel == "significant" else "uniform"
        cfgs.append(CVConfig(name, sel, weighting or default_w, args.coupling, args.combine))
    return cfgs


def cmd_cv(args) -> None:
    mpx, groups = _load(args)
    groups = _select_groups(groups, args)
    configs = _configs(args, mpx)
    _log(f"cv: {len(groups)} groups, configs {[c.name for c in configs]}, "
         f"{args.folds} folds, protocol {args.protocol}")
    t0 = time.time()
    cv = retrieval_cv(mpx, groups, configs, k_folds=args.folds, protocol=args.protocol,
                      r=args.restart, seed=args.seed, min_size=None, max_size=None,
                      lcc_kwargs=_lcc_kwargs(args), progress=True)
    out = Path(args.out)
    cv.folds.to_csv(out / "cv_folds.tsv", sep="\t", index=False)
    cv.held_out.to_csv(out / "cv_held_out.tsv", sep="\t", index=False)
    summary = cv.summary()
    summary.to_csv(out / "cv_summary.tsv", sep="\t", index=False)
    _log(f"cv done ({time.time() - t0:.0f}s)")
    if not args.no_plots:
        from multiome_algo.viz import plot_cv_performance

        plot_cv_performance(cv, out=out / "cv_performance.png")
    print(summary.groupby("config", sort=False)[["auroc_median", "top10", "top100"]]
          .median().to_string())


def cmd_run(args) -> None:
    import yaml

    cfg = yaml.safe_load(Path(args.config).read_text()) or {}
    command = cfg.pop("command", None)
    if command not in COMMANDS or command == "run":
        raise SystemExit(f"config needs 'command: <{'|'.join(c for c in COMMANDS if c != 'run')}>'")
    argv = [command]
    for key, val in cfg.items():
        flag = "--" + key.replace("_", "-")
        if val is True:
            argv.append(flag)
        elif val is False:
            if key in ("normalize_ids",):
                argv.append("--no-" + key.replace("_", "-"))
        elif isinstance(val, list):
            argv += [flag, *map(str, val)]
        elif val is not None:
            argv += [flag, str(val)]
    _log("run: multiome " + " ".join(argv))
    main(argv)


def _safe(s: str) -> str:
    return "".join(c if c.isalnum() or c in "-_." else "_" for c in s)[:80]


COMMANDS = {"modularity": cmd_modularity, "rank": cmd_rank, "cv": cmd_cv, "run": cmd_run}


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="multiome", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="command", required=True)

    p = sub.add_parser("modularity", help="LCC modularity per group x layer")
    _add_input_args(p)
    _add_lcc_args(p)
    p.add_argument("--group", nargs="+", help="only these group ids/labels")
    p.add_argument("--per-group-plots", action="store_true", help="layer plot per group")
    p.add_argument("--no-plots", action="store_true")

    p = sub.add_parser("rank", help="rank candidate genes for one group or seed list")
    _add_input_args(p, groups_required=False)
    _add_lcc_args(p)
    _add_walk_args(p)
    p.add_argument("--group", nargs=1, help="group id/label to use as seeds")
    p.add_argument("--seeds", nargs="+", help="seed genes, or a file with one gene per line")
    p.add_argument("--name", help="name for a --seeds list")
    p.add_argument("--layer-selection", choices=["significant", "all"], default="significant",
                   help="walk significant layers (paper) or all layers uniformly")
    p.add_argument("--modularity", help="reuse a modularity.tsv (paper: BH over all groups)")
    p.add_argument("--highlight", help="file of genes to highlight in the ranking plot")
    p.add_argument("--top", type=int, default=30)
    p.add_argument("--no-plots", action="store_true")

    p = sub.add_parser("cv", help="cross-validated retrieval")
    _add_input_args(p)
    _add_lcc_args(p)
    _add_walk_args(p)
    p.add_argument("--group", nargs="+", help="only these group ids/labels")
    p.add_argument("--folds", type=int, default=10, help="folds per group (paper: 10)")
    p.add_argument("--protocol", choices=["train", "paper"], default="train",
                   help="layer selection from training genes (default) or the full group "
                        "as published")
    p.add_argument("--configs", nargs="+",
                   help="NAME=LAYERS[:WEIGHTING], LAYERS = significant | all | id+id; "
                        "default: informed, all_uniform, ppi")
    p.add_argument("--no-plots", action="store_true")

    p = sub.add_parser("run", help="run a command from a YAML config")
    p.add_argument("config", help="YAML with 'command:' plus option names as keys")
    return ap


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    COMMANDS[args.command](args)


if __name__ == "__main__":
    main()
