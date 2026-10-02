"""Read user data: edge-list folders, layer metadata, and gene sets.

Edge lists
    One file per layer (``.tsv``, ``.csv``, ``.txt``, optionally gzipped); the file stem is
    the layer id. Columns: source, target[, weight]. Delimiter is sniffed. A header row
    is detected automatically (a first row whose genes never occur again is a header) or
    can be forced with ``header=True/False``. Self-loops and missing values are dropped,
    duplicate undirected edges collapsed. A third numeric column is read as edge weight
    only when ``weights=True`` (the paper's layers are all binary).

Layer metadata (optional)
    ``layers.tsv`` in the same folder (or any path passed as ``metadata=``): first column
    is the layer id; optional columns ``tags`` (``;``-separated), ``description``. The
    original repo's ``network_details.tsv`` (columns network/type/subtype/source) is also
    understood.

Gene sets
    * GMT: ``name<TAB>description<TAB>gene1<TAB>gene2...``
    * long table: two columns ``group, gene`` (one row per membership)
    * wide table: one row per group with a ``;``/``,``/``|``-separated gene column
      (e.g. the repo's Orphanet table with columns ``name, all_genes``).
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
import re
from collections.abc import Iterable
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

from multiome_core.schema import GeneGroup, Layer, Multiplex

EDGE_SUFFIXES = (".tsv", ".csv", ".txt", ".tsv.gz", ".csv.gz", ".txt.gz")
METADATA_NAMES = ("layers.tsv", "layers.csv")
DEFAULT_CACHE = Path.home() / ".cache" / "multiome" / "layers"


# --------------------------------------------------------------------------- edge lists


def _open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else open(path)  # noqa: SIM115


def _sniff_delimiter(path: Path) -> str:
    with _open_text(path) as fh:
        lines = [ln for ln in (fh.readline() for _ in range(5)) if ln.strip()]
    for d in ("\t", ",", ";"):  # unambiguous delimiters first: headers may contain spaces
        if lines and all(d in ln for ln in lines):
            return d
    return r"\s+"


def _layer_id(path: Path) -> str:
    name = path.name
    for suf in EDGE_SUFFIXES:
        if name.endswith(suf):
            return name[: -len(suf)]
    return path.stem


def read_edgelist(
    path: str | Path,
    layer_id: str | None = None,
    header: bool | None = None,
    weights: bool = False,
    sep: str | None = None,
    **meta,
) -> Layer:
    """Read one edge-list file into a Layer.

    Args:
        path: edge-list file.
        layer_id: defaults to the file stem.
        header: True/False to force; None to auto-detect.
        weights: read the third column as edge weight.
        sep: delimiter; sniffed when None.
        meta: passed to Layer (tags, description).
    """
    path = Path(path)
    sep = sep or _sniff_delimiter(path)
    ncols = 3 if weights else 2
    df = pd.read_csv(
        path, sep=sep, header=None, usecols=range(ncols), dtype=str,
        engine="c", quoting=csv.QUOTE_NONE,
    )
    df.columns = ["source", "target", "weight"][:ncols]
    for c in ("source", "target"):  # R-written files quote NA as "NA"
        df[c] = df[c].str.strip('"').replace({"NA": None, "": None})

    if header is None:
        first = df.iloc[0]
        rest = df.iloc[1:]
        seen = set(rest["source"].dropna()) | set(rest["target"].dropna())
        header = first["source"] not in seen and first["target"] not in seen
    if header:
        df = df.iloc[1:]

    w = pd.to_numeric(df["weight"], errors="raise") if weights else None
    return Layer.from_edges(
        layer_id or _layer_id(path),
        df["source"],
        df["target"],
        w,
        source=str(path),
        **meta,
    )


def read_layer_metadata(path: str | Path) -> pd.DataFrame:
    """Read layer metadata into columns: layer_id, tags (tuple), description."""
    path = Path(path)
    df = pd.read_csv(path, sep=None, engine="python", dtype=str).fillna("")
    cols = {c.lower(): c for c in df.columns}
    out = pd.DataFrame({"layer_id": df.iloc[:, 0]})
    if "tags" in cols:
        out["tags"] = df[cols["tags"]].map(lambda s: tuple(t for t in s.split(";") if t))
    elif "type" in cols:  # original network_details.tsv: type, else main_type
        main = df[cols["main_type"]] if "main_type" in cols else pd.Series("", index=df.index)
        out["tags"] = [(t or m,) if (t or m) else () for t, m in zip(df[cols["type"]], main, strict=True)]
    else:
        out["tags"] = [()] * len(df)
    desc_col = cols.get("description") or cols.get("subtype")
    out["description"] = df[desc_col] if desc_col else ""
    return out


def _cache_key(path: Path, header, weights) -> str:
    st = path.stat()
    raw = f"{path.resolve()}|{st.st_size}|{st.st_mtime_ns}|{header}|{weights}|v2"
    return hashlib.sha1(raw.encode()).hexdigest()[:16]


def _load_cached(path: Path, cache_dir: Path | None, header, weights, **meta) -> Layer:
    if cache_dir is None:
        return read_edgelist(path, header=header, weights=weights, **meta)
    f = cache_dir / f"{_layer_id(path)}-{_cache_key(path, header, weights)}.npz"
    if f.exists():
        z = np.load(f, allow_pickle=True)
        adj = sp.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"]))
        return Layer(
            id=_layer_id(path), genes=z["genes"], adj=adj, weighted=bool(z["weighted"]),
            source=str(path), **meta,
        )
    layer = read_edgelist(path, header=header, weights=weights, **meta)
    cache_dir.mkdir(parents=True, exist_ok=True)
    np.savez(
        f, data=layer.adj.data, indices=layer.adj.indices, indptr=layer.adj.indptr,
        shape=np.array(layer.adj.shape), genes=layer.genes, weighted=layer.weighted,
    )
    return layer


def read_multiplex(
    folder: str | Path,
    layer_ids: Iterable[str] | None = None,
    metadata: str | Path | None = None,
    header: bool | None = None,
    weights: bool = False,
    cache_dir: str | Path | None = DEFAULT_CACHE,
    name: str = "",
) -> Multiplex:
    """Read a folder of edge lists (one file per layer) into a Multiplex.

    Args:
        folder: directory of edge-list files.
        layer_ids: only load these layers (file stems).
        metadata: layer metadata table; defaults to ``layers.tsv``/``layers.csv`` in folder.
        header, weights: see `read_edgelist`.
        cache_dir: parsed layers are cached here as sparse matrices (None disables).
        name: name for the multiplex.
    """
    folder = Path(folder)
    files = {
        _layer_id(p): p
        for p in sorted(folder.iterdir())
        if p.name.endswith(EDGE_SUFFIXES) and p.name not in METADATA_NAMES
    }
    if layer_ids is not None:
        layer_ids = list(layer_ids)
        missing = [lid for lid in layer_ids if lid not in files]
        if missing:
            raise FileNotFoundError(f"No edge list for layers {missing} in {folder}")
        files = {lid: files[lid] for lid in layer_ids}
    if not files:
        raise FileNotFoundError(f"No edge-list files ({', '.join(EDGE_SUFFIXES)}) in {folder}")

    meta: dict[str, dict] = {}
    if metadata is None:
        metadata = next((folder / m for m in METADATA_NAMES if (folder / m).exists()), None)
    if metadata is not None:
        for row in read_layer_metadata(metadata).itertuples(index=False):
            meta[row.layer_id] = {"tags": row.tags, "description": row.description}

    cache = Path(cache_dir) if cache_dir is not None else None
    mpx = Multiplex(name=name or folder.name)
    for lid, path in files.items():
        mpx.add(_load_cached(path, cache, header, weights, **meta.get(lid, {})))
    return mpx


def multiplex_from_dataframes(
    edges: dict[str, pd.DataFrame],
    source: str = "source",
    target: str = "target",
    weight: str | None = None,
    name: str = "",
) -> Multiplex:
    """Build a Multiplex from in-memory edge tables ({layer_id: DataFrame})."""
    mpx = Multiplex(name=name)
    for lid, df in edges.items():
        w = df[weight] if weight and weight in df else None
        mpx.add(Layer.from_edges(lid, df[source], df[target], w))
    return mpx


# --------------------------------------------------------------------------- gene sets

_GENE_SPLIT = re.compile(r"[;,|]\s*")


def _split_genes(s: str) -> list[str]:
    return [g.strip() for g in _GENE_SPLIT.split(s) if g.strip()]


def _read_gmt(text: str, source: str) -> list[GeneGroup]:
    groups = []
    for line in text.splitlines():
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 3 or not parts[0]:
            continue
        genes = [g.strip() for g in parts[2:] if g.strip()]
        groups.append(
            GeneGroup(id=parts[0], genes=genes, label=parts[0], source=source,
                      metadata={"description": parts[1]})
        )
    return groups


def read_gene_groups(
    path: str | Path,
    fmt: str = "auto",
    id_col: str | None = None,
    label_col: str | None = None,
    genes_col: str | None = None,
) -> list[GeneGroup]:
    """Read gene groups from a GMT, long (group, gene) or wide (one row per group) table.

    Args:
        path: input file.
        fmt: "gmt", "long", "wide" or "auto" (by extension, then by content).
        id_col: column holding a stable group id (wide/long). Auto: an ``id``-like column
            (e.g. ``orphaID``), else the label.
        label_col: column with the human-readable name. Auto: ``name``/``label``.
        genes_col: wide format - column holding the separated gene list.
            Auto: ``all_genes``/``genes``.
    """
    path = Path(path)
    with _open_text(path) as fh:
        text = fh.read()
    source = path.name

    if fmt == "auto":
        fmt = "gmt" if path.name.endswith((".gmt", ".gmt.gz")) else ""
    if fmt == "gmt":
        return _read_gmt(text, source)

    df = pd.read_csv(io.StringIO(text), sep=None, engine="python", dtype=str)
    cols = {c.lower(): c for c in df.columns}

    def pick(explicit, *candidates):
        if explicit:
            return explicit
        return next((cols[c] for c in candidates if c in cols), None)

    genes_c = pick(genes_col, "all_genes", "genes", "gene_list", "members")
    if fmt == "" and genes_c is None and df.shape[1] == 2:
        fmt = "long"
    if fmt in ("", "wide") and genes_c is not None:
        label_c = pick(label_col, "name", "label", "group", "term") or df.columns[0]
        id_c = pick(id_col, "id", "group_id", "orphaid", "term_id") or label_c
        df = df.dropna(subset=[genes_c])
        groups = []
        seen: dict[str, int] = {}
        for r in df.to_dict("records"):
            gid = str(r[id_c])
            if gid in seen:  # keep ids unique; duplicated labels exist in real tables
                seen[gid] += 1
                gid = f"{gid}#{seen[gid]}"
            else:
                seen[gid] = 0
            groups.append(
                GeneGroup(id=gid, genes=_split_genes(r[genes_c]), label=str(r[label_c]),
                          source=source)
            )
        return groups
    if fmt == "long":
        grp_c = pick(id_col, "group", "group_id", "set", "term") or df.columns[0]
        gene_c = pick(genes_col, "gene", "symbol", "gene_symbol") or df.columns[1]
        df = df.dropna(subset=[grp_c, gene_c])
        return [
            GeneGroup(id=str(gid), genes=sub[gene_c].str.strip(), source=source)
            for gid, sub in df.groupby(grp_c, sort=False)
        ]
    raise ValueError(
        f"Could not infer gene-set format of {path}; pass fmt='gmt'|'long'|'wide' "
        "and the relevant column names."
    )


def write_gmt(groups: Iterable[GeneGroup], path: str | Path) -> None:
    with open(path, "w") as fh:
        for g in groups:
            fh.write("\t".join([g.id, g.label, *sorted(g.genes)]) + "\n")


# --------------------------------------------------------------------------- coverage


def coverage_report(multiplex: Multiplex, groups: Iterable[GeneGroup]) -> pd.DataFrame:
    """Per group: how many genes are found in the multiplex and in each layer.

    Columns: group_id, n_genes, n_in_universe, n_missing, missing (;-joined),
    and one ``n_in:<layer>`` column per layer.
    """
    universe = multiplex.gene_set()
    rows = []
    for g in groups:
        missing = sorted(g.genes - universe)
        row = {
            "group_id": g.id,
            "label": g.label,
            "n_genes": len(g),
            "n_in_universe": len(g) - len(missing),
            "n_missing": len(missing),
            "missing": ";".join(missing),
        }
        for layer in multiplex:
            row[f"n_in:{layer.id}"] = len(g.genes & layer.gene_set)
        rows.append(row)
    return pd.DataFrame(rows)
