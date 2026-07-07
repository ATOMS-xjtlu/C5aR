#!/usr/bin/env python3
"""NetworkX key-residue route analysis for C5aR1 AlloViz outputs.

This script follows the AlloViz NetworkX tutorial pattern on the already
exported formal networks. It reconstructs NetworkX graphs from the filtered
edge tables, computes shortest paths between curated C5aR1 key residues, and
summarizes repeat-level path frequency and consensus routes.
"""

from __future__ import annotations

import itertools
import math
import os
import re
from collections import Counter, defaultdict
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-alloviz-c5ar")

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RUNS_DIR = ROOT / "runs" / "formal"
TABLES_DIR = ROOT / "tables"
FIGURES_DIR = ROOT / "figures"
LOGS_DIR = ROOT / "logs"

MAPPING_CSV = TABLES_DIR / "residue_mapping_alloviz.csv"
OUT_PATHS = TABLES_DIR / "networkx_key_shortest_paths_long.csv"
OUT_PATH_FREQ = TABLES_DIR / "networkx_key_exact_path_frequency.csv"
OUT_NODE_FREQ = TABLES_DIR / "networkx_key_path_node_frequency.csv"
OUT_EDGE_FREQ = TABLES_DIR / "networkx_key_path_edge_frequency.csv"
OUT_GLOBAL_NODE_FREQ = TABLES_DIR / "networkx_global_path_node_frequency.csv"
OUT_GLOBAL_EDGE_FREQ = TABLES_DIR / "networkx_global_path_edge_frequency.csv"
OUT_CONSENSUS = TABLES_DIR / "networkx_consensus_routes.csv"
OUT_FOCUS_COMPARE = TABLES_DIR / "networkx_focus_pair_distance_comparison.csv"
OUT_REPORT = ROOT / "networkx_key_route_report.md"
OUT_HEATMAP = FIGURES_DIR / "networkx_key_route_distance_delta_heatmap.png"
OUT_EDGE_BAR = FIGURES_DIR / "networkx_global_consensus_edge_support.png"
OUT_STATUS = LOGS_DIR / "networkx_key_routes_status.json"

SYSTEMS = ["apo", "BM213", "C5apep"]
REPLICATES = ["rep1", "rep2", "rep3"]
PKGS = ["GetContacts", "correlationplus_CA_Pear", "pytraj_CA"]
FILTERING = "Spatially_distant"
EDGE_EPS = 1e-9

EDGE_NODE_RE = re.compile(r"([A-Z]{3}):(-?\d+)")

KEY_RESIDUES_ORDER = [
    "I91",
    "W102",
    "I116",
    "M120",
    "R134",
    "S171",
    "Y222",
    "F251",
    "W255",
    "Y258",
    "Y290",
    "N292",
    "N296",
    "Y300",
]

FOCUS_PAIRS = [
    ("W102", "M120"),
    ("I116", "M120"),
    ("I116", "Y258"),
    ("S171", "Y258"),
    ("M120", "W255"),
    ("M120", "Y258"),
    ("R134", "Y222"),
    ("Y222", "F251"),
    ("F251", "W255"),
    ("W255", "Y258"),
    ("W255", "Y290"),
    ("Y258", "Y300"),
    ("M120", "Y300"),
    ("N296", "Y300"),
]


def parse_edge_node(value: str) -> tuple[str, int]:
    match = EDGE_NODE_RE.search(str(value))
    if not match:
        raise ValueError(f"Cannot parse AlloViz edge node label: {value!r}")
    return match.group(1), int(match.group(2))


def markdown_table(rows: list[dict], columns: list[str], max_rows: int | None = None) -> str:
    if max_rows is not None:
        rows = rows[:max_rows]
    if not rows:
        return "_No rows._"
    header = "| " + " | ".join(columns) + " |"
    sep = "| " + " | ".join(["---"] * len(columns)) + " |"
    body = []
    for row in rows:
        vals = []
        for col in columns:
            value = row.get(col, "")
            if isinstance(value, float):
                if math.isnan(value):
                    value = ""
                elif value.is_integer() and (
                    col.endswith("count")
                    or col.endswith("counts")
                    or col in {"ok", "missing_node", "no_path", "graph_error", "total_pairs_x_repeats", "support_repeats", "successful_repeats", "modal_path_support_repeats", "path_occurrences", "unique_pairs", "focus_pair_occurrences", "intermediate_occurrences"}
                ):
                    value = f"{int(value)}"
                else:
                    value = f"{value:.3f}"
            vals.append(str(value).replace("|", "/"))
        body.append("| " + " | ".join(vals) + " |")
    return "\n".join([header, sep, *body])


def load_mapping() -> tuple[dict, dict, dict]:
    mapping = pd.read_csv(MAPPING_CSV)
    mapping["biological_residue"] = mapping["biological_residue"].fillna("")
    mapping["gpcrdb"] = mapping["gpcrdb"].fillna("")
    mapping["panel_context"] = mapping["panel_context"].fillna("")

    raw_to_canonical: dict[str, dict[int, int]] = {}
    canonical_info: dict[int, dict] = {}
    key_to_canonical: dict[str, int] = {}

    for system, sdf in mapping.groupby("system"):
        raw_to_canonical[system] = {
            int(row.raw_resid): int(row.canonical_resid) for row in sdf.itertuples()
        }

    for row in mapping.sort_values(["system", "canonical_resid"]).itertuples():
        canonical = int(row.canonical_resid)
        if canonical not in canonical_info:
            inferred_number = canonical + 33
            canonical_info[canonical] = {
                "canonical_resid": canonical,
                "canonical_label": row.canonical_label,
                "aa": row.aa,
                "display_label": f"{row.aa}{inferred_number}",
                "biological_residue": "",
                "gpcrdb": "",
                "panel_context": "",
                "is_key_residue": False,
            }
        if row.biological_residue:
            canonical_info[canonical].update(
                {
                    "display_label": row.biological_residue,
                    "biological_residue": row.biological_residue,
                    "gpcrdb": row.gpcrdb,
                    "panel_context": row.panel_context,
                    "is_key_residue": bool(row.is_key_residue),
                }
            )
            key_to_canonical[row.biological_residue] = canonical

    return raw_to_canonical, canonical_info, key_to_canonical


def node_label(canonical: int, canonical_info: dict[int, dict]) -> str:
    return canonical_info.get(canonical, {}).get("display_label", f"X{canonical + 33}")


def node_context(canonical: int, canonical_info: dict[int, dict]) -> str:
    return canonical_info.get(canonical, {}).get("panel_context", "")


def edge_label(a: int, b: int, canonical_info: dict[int, dict]) -> str:
    aa, bb = sorted((a, b))
    return f"{node_label(aa, canonical_info)}-{node_label(bb, canonical_info)}"


def path_to_labels(path: list[int], canonical_info: dict[int, dict]) -> str:
    return " -> ".join(node_label(n, canonical_info) for n in path)


def path_to_canonicals(path: list[int]) -> str:
    return " -> ".join(str(n) for n in path)


def path_edge_labels(path: list[int], canonical_info: dict[int, dict]) -> str:
    return "; ".join(edge_label(a, b, canonical_info) for a, b in zip(path[:-1], path[1:]))


def edge_csv_path(system: str, rep: str, pkg: str) -> Path:
    return (
        RUNS_DIR
        / system
        / rep
        / f"{pkg}_serial_c1_downsampled"
        / "exports"
        / f"{pkg}__{FILTERING}__edges.csv"
    )


def build_graph(system: str, rep: str, pkg: str, raw_to_canonical: dict) -> tuple[nx.Graph, dict]:
    path = edge_csv_path(system, rep, pkg)
    if not path.exists():
        raise FileNotFoundError(path)

    df = pd.read_csv(path)
    graph = nx.Graph()
    skipped = Counter()

    for row in df.itertuples(index=False):
        left = getattr(row, "_0")
        right = getattr(row, "_1")
        try:
            _resname_a, raw_a = parse_edge_node(left)
            _resname_b, raw_b = parse_edge_node(right)
        except ValueError:
            skipped["parse_error"] += 1
            continue

        if raw_a not in raw_to_canonical[system] or raw_b not in raw_to_canonical[system]:
            skipped["unmapped_raw_resid"] += 1
            continue

        ca = raw_to_canonical[system][raw_a]
        cb = raw_to_canonical[system][raw_b]
        if ca == cb:
            skipped["self_edge"] += 1
            continue

        weight = abs(float(row.weight))
        if not np.isfinite(weight):
            skipped["nonfinite_weight"] += 1
            continue

        distance = -math.log(weight + EDGE_EPS) + EDGE_EPS
        if distance < 0:
            distance = 0.0

        if graph.has_edge(ca, cb):
            if weight <= graph[ca][cb]["graph_weight"]:
                continue

        graph.add_edge(
            ca,
            cb,
            graph_weight=weight,
            graph_distance=distance,
            raw_weight=weight,
            source_csv=str(path),
        )

    meta = {
        "edge_csv": str(path),
        "edges_in_csv": int(len(df)),
        "nodes_in_graph": int(graph.number_of_nodes()),
        "edges_in_graph": int(graph.number_of_edges()),
        "skipped": dict(skipped),
    }
    return graph, meta


def get_key_pairs(key_to_canonical: dict[str, int]) -> tuple[list[tuple[int, int]], set[tuple[int, int]]]:
    ordered = [key_to_canonical[k] for k in KEY_RESIDUES_ORDER if k in key_to_canonical]
    pairs = [tuple(sorted(pair)) for pair in itertools.combinations(ordered, 2)]
    focus = set()
    for left, right in FOCUS_PAIRS:
        if left in key_to_canonical and right in key_to_canonical:
            focus.add(tuple(sorted((key_to_canonical[left], key_to_canonical[right]))))
    return pairs, focus


def compute_shortest_paths() -> tuple[pd.DataFrame, list[dict]]:
    raw_to_canonical, canonical_info, key_to_canonical = load_mapping()
    key_pairs, focus_pairs = get_key_pairs(key_to_canonical)

    records: list[dict] = []
    graph_status: list[dict] = []

    for system in SYSTEMS:
        for rep in REPLICATES:
            for pkg in PKGS:
                try:
                    graph, meta = build_graph(system, rep, pkg, raw_to_canonical)
                    graph_ok = True
                except Exception as exc:
                    graph = nx.Graph()
                    meta = {
                        "edge_csv": str(edge_csv_path(system, rep, pkg)),
                        "edges_in_csv": 0,
                        "nodes_in_graph": 0,
                        "edges_in_graph": 0,
                        "skipped": {},
                        "error": repr(exc),
                    }
                    graph_ok = False

                graph_status.append({"system": system, "replicate": rep, "pkg": pkg, **meta})

                for src, dst in key_pairs:
                    source_label = node_label(src, canonical_info)
                    target_label = node_label(dst, canonical_info)
                    pair_label = f"{source_label}-{target_label}"
                    base = {
                        "system": system,
                        "replicate": rep,
                        "pkg": pkg,
                        "filtering": FILTERING,
                        "source_canonical_resid": src,
                        "target_canonical_resid": dst,
                        "source_residue": source_label,
                        "target_residue": target_label,
                        "pair_label": pair_label,
                        "is_focus_pair": (tuple(sorted((src, dst))) in focus_pairs),
                    }

                    if not graph_ok:
                        records.append({**base, "status": "graph_error"})
                        continue
                    if src not in graph or dst not in graph:
                        missing = []
                        if src not in graph:
                            missing.append(source_label)
                        if dst not in graph:
                            missing.append(target_label)
                        records.append(
                            {**base, "status": "missing_node", "missing_nodes": ";".join(missing)}
                        )
                        continue
                    if not nx.has_path(graph, src, dst):
                        records.append({**base, "status": "no_path"})
                        continue

                    path = nx.shortest_path(graph, src, dst, weight="graph_distance")
                    distance = nx.shortest_path_length(graph, src, dst, weight="graph_distance")
                    weights = [
                        graph[a][b]["graph_weight"] for a, b in zip(path[:-1], path[1:])
                    ]
                    records.append(
                        {
                            **base,
                            "status": "ok",
                            "path_length_edges": len(path) - 1,
                            "path_distance": float(distance),
                            "path_mean_graph_weight": float(np.mean(weights)) if weights else np.nan,
                            "path_min_graph_weight": float(np.min(weights)) if weights else np.nan,
                            "path_canonical_resids": path_to_canonicals(path),
                            "path_residues": path_to_labels(path, canonical_info),
                            "path_edges": path_edge_labels(path, canonical_info),
                            "intermediate_residues": path_to_labels(path[1:-1], canonical_info)
                            if len(path) > 2
                            else "",
                        }
                    )

    return pd.DataFrame.from_records(records), graph_status


def explode_path_nodes(paths: pd.DataFrame, canonical_info: dict[int, dict]) -> pd.DataFrame:
    rows = []
    for rec in paths[paths["status"].eq("ok")].itertuples(index=False):
        nodes = [int(x) for x in str(rec.path_canonical_resids).split(" -> ")]
        for order, node in enumerate(nodes):
            rows.append(
                {
                    "system": rec.system,
                    "replicate": rec.replicate,
                    "pkg": rec.pkg,
                    "source_residue": rec.source_residue,
                    "target_residue": rec.target_residue,
                    "pair_label": rec.pair_label,
                    "is_focus_pair": bool(rec.is_focus_pair),
                    "node_order": order,
                    "node_canonical_resid": node,
                    "node_residue": node_label(node, canonical_info),
                    "node_context": node_context(node, canonical_info),
                    "is_endpoint": order == 0 or order == len(nodes) - 1,
                    "is_intermediate": 0 < order < len(nodes) - 1,
                }
            )
    return pd.DataFrame.from_records(rows)


def explode_path_edges(paths: pd.DataFrame, canonical_info: dict[int, dict]) -> pd.DataFrame:
    rows = []
    for rec in paths[paths["status"].eq("ok")].itertuples(index=False):
        nodes = [int(x) for x in str(rec.path_canonical_resids).split(" -> ")]
        for order, (a, b) in enumerate(zip(nodes[:-1], nodes[1:])):
            aa, bb = sorted((a, b))
            rows.append(
                {
                    "system": rec.system,
                    "replicate": rec.replicate,
                    "pkg": rec.pkg,
                    "source_residue": rec.source_residue,
                    "target_residue": rec.target_residue,
                    "pair_label": rec.pair_label,
                    "is_focus_pair": bool(rec.is_focus_pair),
                    "edge_order": order,
                    "edge_canonical_a": aa,
                    "edge_canonical_b": bb,
                    "edge_label": edge_label(aa, bb, canonical_info),
                }
            )
    return pd.DataFrame.from_records(rows)


def summarize(paths: pd.DataFrame, canonical_info: dict[int, dict]) -> dict[str, pd.DataFrame]:
    ok = paths[paths["status"].eq("ok")].copy()

    exact_path_freq = (
        ok.groupby(
            [
                "system",
                "pkg",
                "source_residue",
                "target_residue",
                "pair_label",
                "is_focus_pair",
                "path_residues",
                "path_canonical_resids",
            ],
            dropna=False,
        )
        .agg(
            support_repeats=("replicate", "nunique"),
            replicates=("replicate", lambda s: ",".join(sorted(set(s)))),
            mean_path_distance=("path_distance", "mean"),
            sd_path_distance=("path_distance", "std"),
            mean_path_length_edges=("path_length_edges", "mean"),
        )
        .reset_index()
        .sort_values(
            [
                "system",
                "pkg",
                "is_focus_pair",
                "support_repeats",
                "mean_path_distance",
            ],
            ascending=[True, True, False, False, True],
        )
    )

    node_occ = explode_path_nodes(paths, canonical_info)
    edge_occ = explode_path_edges(paths, canonical_info)

    node_freq = (
        node_occ.groupby(
            [
                "system",
                "pkg",
                "source_residue",
                "target_residue",
                "pair_label",
                "is_focus_pair",
                "node_canonical_resid",
                "node_residue",
                "node_context",
                "is_endpoint",
                "is_intermediate",
            ],
            dropna=False,
        )
        .agg(
            support_repeats=("replicate", "nunique"),
            replicates=("replicate", lambda s: ",".join(sorted(set(s)))),
            occurrences=("replicate", "size"),
        )
        .reset_index()
        .sort_values(["system", "pkg", "pair_label", "support_repeats"], ascending=[True, True, True, False])
    )

    edge_freq = (
        edge_occ.groupby(
            [
                "system",
                "pkg",
                "source_residue",
                "target_residue",
                "pair_label",
                "is_focus_pair",
                "edge_canonical_a",
                "edge_canonical_b",
                "edge_label",
            ],
            dropna=False,
        )
        .agg(
            support_repeats=("replicate", "nunique"),
            replicates=("replicate", lambda s: ",".join(sorted(set(s)))),
            occurrences=("replicate", "size"),
        )
        .reset_index()
        .sort_values(["system", "pkg", "pair_label", "support_repeats"], ascending=[True, True, True, False])
    )

    global_node_freq = (
        node_occ.groupby(["system", "pkg", "node_canonical_resid", "node_residue", "node_context"], dropna=False)
        .agg(
            path_occurrences=("pair_label", "size"),
            unique_pairs=("pair_label", "nunique"),
            unique_pair_repeats=("replicate", "size"),
            intermediate_occurrences=("is_intermediate", "sum"),
            endpoint_occurrences=("is_endpoint", "sum"),
            focus_pair_occurrences=("is_focus_pair", "sum"),
        )
        .reset_index()
        .sort_values(["system", "pkg", "intermediate_occurrences", "path_occurrences"], ascending=[True, True, False, False])
    )

    global_edge_freq = (
        edge_occ.groupby(["system", "pkg", "edge_canonical_a", "edge_canonical_b", "edge_label"], dropna=False)
        .agg(
            path_occurrences=("pair_label", "size"),
            unique_pairs=("pair_label", "nunique"),
            support_replicates=("replicate", "nunique"),
            focus_pair_occurrences=("is_focus_pair", "sum"),
        )
        .reset_index()
        .sort_values(["system", "pkg", "path_occurrences", "unique_pairs"], ascending=[True, True, False, False])
    )

    consensus_rows = []
    for keys, group in ok.groupby(
        ["system", "pkg", "source_residue", "target_residue", "pair_label", "is_focus_pair"],
        dropna=False,
    ):
        system, pkg, source, target, pair_label, is_focus_pair = keys
        exact = (
            group.groupby(["path_residues", "path_canonical_resids"], dropna=False)
            .agg(
                support_repeats=("replicate", "nunique"),
                replicates=("replicate", lambda s: ",".join(sorted(set(s)))),
                mean_path_distance=("path_distance", "mean"),
                mean_path_length_edges=("path_length_edges", "mean"),
            )
            .reset_index()
            .sort_values(["support_repeats", "mean_path_distance"], ascending=[False, True])
        )
        modal = exact.iloc[0]

        pair_node_freq = node_freq[
            (node_freq["system"].eq(system))
            & (node_freq["pkg"].eq(pkg))
            & (node_freq["pair_label"].eq(pair_label))
        ]
        pair_edge_freq = edge_freq[
            (edge_freq["system"].eq(system))
            & (edge_freq["pkg"].eq(pkg))
            & (edge_freq["pair_label"].eq(pair_label))
        ]
        consensus_nodes = pair_node_freq[pair_node_freq["support_repeats"].ge(2)]
        consensus_intermediates = consensus_nodes[consensus_nodes["is_intermediate"].astype(bool)]
        consensus_edges = pair_edge_freq[pair_edge_freq["support_repeats"].ge(2)]

        consensus_rows.append(
            {
                "system": system,
                "pkg": pkg,
                "source_residue": source,
                "target_residue": target,
                "pair_label": pair_label,
                "is_focus_pair": bool(is_focus_pair),
                "successful_repeats": int(group["replicate"].nunique()),
                "mean_path_distance": float(group["path_distance"].mean()),
                "sd_path_distance": float(group["path_distance"].std())
                if group["path_distance"].notna().sum() > 1
                else np.nan,
                "mean_path_length_edges": float(group["path_length_edges"].mean()),
                "modal_path_support_repeats": int(modal.support_repeats),
                "modal_path_replicates": modal.replicates,
                "modal_path_distance_mean": float(modal.mean_path_distance),
                "modal_path_residues": modal.path_residues,
                "consensus_nodes_2of3": "; ".join(consensus_nodes["node_residue"].tolist()),
                "consensus_intermediate_nodes_2of3": "; ".join(
                    consensus_intermediates["node_residue"].tolist()
                ),
                "consensus_edges_2of3": "; ".join(consensus_edges["edge_label"].tolist()),
            }
        )

    consensus = pd.DataFrame.from_records(consensus_rows).sort_values(
        ["system", "pkg", "is_focus_pair", "modal_path_support_repeats", "mean_path_distance"],
        ascending=[True, True, False, False, True],
    )

    focus = consensus[consensus["is_focus_pair"].astype(bool)].copy()
    focus_compare = (
        focus.pivot_table(
            index=["pkg", "pair_label"],
            columns="system",
            values="mean_path_distance",
            aggfunc="mean",
        )
        .reset_index()
        .rename_axis(None, axis=1)
    )
    for col in SYSTEMS:
        if col not in focus_compare:
            focus_compare[col] = np.nan
    focus_compare["delta_BM213_minus_apo"] = focus_compare["BM213"] - focus_compare["apo"]
    focus_compare["delta_C5apep_minus_apo"] = focus_compare["C5apep"] - focus_compare["apo"]
    focus_compare["delta_BM213_minus_C5apep"] = focus_compare["BM213"] - focus_compare["C5apep"]
    focus_compare = focus_compare.sort_values(["pkg", "delta_BM213_minus_C5apep"], na_position="last")

    return {
        "exact_path_freq": exact_path_freq,
        "node_freq": node_freq,
        "edge_freq": edge_freq,
        "global_node_freq": global_node_freq,
        "global_edge_freq": global_edge_freq,
        "consensus": consensus,
        "focus_compare": focus_compare,
    }


def plot_focus_heatmap(focus_compare: pd.DataFrame) -> None:
    pkgs = [pkg for pkg in PKGS if pkg in set(focus_compare["pkg"])]
    if not pkgs:
        return
    pairs = [f"{a}-{b}" for a, b in FOCUS_PAIRS]
    fig, axes = plt.subplots(1, len(pkgs), figsize=(5.9 * len(pkgs), 9.0), sharey=True)
    if len(pkgs) == 1:
        axes = [axes]
    cmap = plt.get_cmap("coolwarm")

    for ax, pkg in zip(axes, pkgs):
        sub = focus_compare[focus_compare["pkg"].eq(pkg)].set_index("pair_label")
        values = np.array([[sub.loc[p, "delta_BM213_minus_C5apep"] if p in sub.index else np.nan] for p in pairs])
        finite = values[np.isfinite(values)]
        vmax = max(0.1, float(np.nanmax(np.abs(finite))) if finite.size else 1.0)
        im = ax.imshow(values, aspect="auto", cmap=cmap, vmin=-vmax, vmax=vmax)
        ax.set_title(pkg, fontsize=16, pad=10)
        ax.set_xticks([0])
        ax.set_xticklabels(["BM213 - C5apep"], fontsize=18)
        ax.set_yticks(range(len(pairs)))
        ax.set_yticklabels(pairs, fontsize=18)
        for y, val in enumerate(values[:, 0]):
            label = "NA" if np.isnan(val) else f"{val:.2f}"
            ax.text(0, y, label, ha="center", va="center", fontsize=12)
        ax.tick_params(axis="x", labelrotation=0, pad=8)
        cbar = fig.colorbar(im, ax=ax, fraction=0.08, pad=0.02)
        cbar.set_label("Delta distance", fontsize=13)
        cbar.ax.tick_params(labelsize=11)

    fig.suptitle("NetworkX shortest-path distance deltas for C5aR1 focus pairs", fontsize=18, y=0.98)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_HEATMAP, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_edge_support(global_edge_freq: pd.DataFrame) -> None:
    rows = []
    for system in ["BM213", "C5apep", "apo"]:
        for pkg in PKGS:
            sub = global_edge_freq[
                global_edge_freq["system"].eq(system) & global_edge_freq["pkg"].eq(pkg)
            ].head(8)
            for rec in sub.itertuples(index=False):
                rows.append(
                    {
                        "system_pkg": f"{system}\n{pkg}",
                        "edge_label": rec.edge_label,
                        "path_occurrences": rec.path_occurrences,
                    }
                )
    if not rows:
        return
    df = pd.DataFrame(rows)
    # Show the globally most repeated edges in a compact horizontal chart.
    top_edges = (
        df.groupby("edge_label")["path_occurrences"]
        .sum()
        .sort_values(ascending=False)
        .head(20)
        .index.tolist()
    )
    plot_df = df[df["edge_label"].isin(top_edges)].copy()
    edge_totals = plot_df.groupby("edge_label")["path_occurrences"].sum().sort_values()

    fig, ax = plt.subplots(figsize=(10, 7))
    y = np.arange(len(edge_totals))
    ax.barh(y, edge_totals.values, color="#5A8F8C")
    ax.set_yticks(y)
    ax.set_yticklabels(edge_totals.index)
    ax.set_xlabel("Path occurrences across key-pair shortest paths")
    ax.set_title("Most frequent shortest-path edges across C5aR1 key-residue routes")
    fig.tight_layout()
    fig.savefig(OUT_EDGE_BAR, dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_report(paths: pd.DataFrame, summaries: dict[str, pd.DataFrame], graph_status: list[dict]) -> None:
    consensus = summaries["consensus"]
    focus_compare = summaries["focus_compare"]
    global_node = summaries["global_node_freq"]
    global_edge = summaries["global_edge_freq"]
    exact = summaries["exact_path_freq"]

    status = pd.DataFrame(graph_status)
    coverage = (
        paths.groupby(["system", "pkg", "status"])
        .size()
        .reset_index(name="count")
        .pivot_table(index=["system", "pkg"], columns="status", values="count", fill_value=0)
        .reset_index()
    )
    for col in ["ok", "missing_node", "no_path", "graph_error"]:
        if col not in coverage:
            coverage[col] = 0
    coverage["total_pairs_x_repeats"] = coverage[["ok", "missing_node", "no_path", "graph_error"]].sum(axis=1)
    coverage["ok_fraction"] = coverage["ok"] / coverage["total_pairs_x_repeats"]
    coverage_rows = coverage.sort_values(["system", "pkg"]).to_dict("records")

    bm_shorter = (
        focus_compare[focus_compare["delta_BM213_minus_C5apep"].notna()]
        .sort_values("delta_BM213_minus_C5apep")
        .head(12)
        .to_dict("records")
    )
    c5_shorter = (
        focus_compare[focus_compare["delta_BM213_minus_C5apep"].notna()]
        .sort_values("delta_BM213_minus_C5apep", ascending=False)
        .head(12)
        .to_dict("records")
    )

    repeated_focus = (
        exact[(exact["is_focus_pair"].astype(bool)) & (exact["support_repeats"].ge(2))]
        .sort_values(["support_repeats", "mean_path_distance"], ascending=[False, True])
        .head(20)
        .to_dict("records")
    )

    top_nodes = (
        global_node[global_node["intermediate_occurrences"].gt(0)]
        .sort_values(["pkg", "system", "intermediate_occurrences"], ascending=[True, True, False])
        .groupby(["system", "pkg"])
        .head(5)
        .to_dict("records")
    )

    top_edges = (
        global_edge.sort_values(["system", "pkg", "path_occurrences"], ascending=[True, True, False])
        .groupby(["system", "pkg"])
        .head(5)
        .to_dict("records")
    )

    consensus_focus = (
        consensus[consensus["is_focus_pair"].astype(bool)]
        .sort_values(["system", "pkg", "modal_path_support_repeats", "mean_path_distance"], ascending=[True, True, False, True])
        .groupby(["system", "pkg"])
        .head(4)
        .to_dict("records")
    )

    lines = [
        "# C5aR1 NetworkX key-residue route analysis",
        "",
        "## ",
        "",
        " AlloViz  AlloViz  NetworkX graph",
        " key residues  shortest path repeat  consensus route",
        "",
        "- `runs/formal/*/*/*_serial_c1_downsampled/exports/*__Spatially_distant__edges.csv`",
        "- GetContactscorrelationplus_CA_Pearpytraj_CA 3  repeat ",
        "-  `canonical_resid`  C5aR1 normal numberingkey residues  mapping  `biological_residue` key  offset  `AA(canonical_resid+33)`",
        "-  AlloViz graph `graph_distance = -log(abs(weight) + 1e-9) + 1e-9` aggregate network connection ",
        "- shortest path ",
        "",
        "## ",
        "",
        "|  |  |",
        "|---|---|",
        f"| `{OUT_PATHS.relative_to(ROOT)}` |  system/replicate/pkg/key-pair  shortest path  |",
        f"| `{OUT_PATH_FREQ.relative_to(ROOT)}` | exact path  3  repeat  |",
        f"| `{OUT_NODE_FREQ.relative_to(ROOT)}` |  key pair  path node repeat frequency |",
        f"| `{OUT_EDGE_FREQ.relative_to(ROOT)}` |  key pair  path edge repeat frequency |",
        f"| `{OUT_GLOBAL_NODE_FREQ.relative_to(ROOT)}` |  key-pair paths  |",
        f"| `{OUT_GLOBAL_EDGE_FREQ.relative_to(ROOT)}` |  key-pair paths  |",
        f"| `{OUT_CONSENSUS.relative_to(ROOT)}` |  key pair  modal path  2/3 consensus nodes/edges |",
        f"| `{OUT_FOCUS_COMPARE.relative_to(ROOT)}` |  pair  BM213/C5apep/apo path distance  |",
        f"| `{OUT_HEATMAP.relative_to(ROOT)}` | BM213 - C5apep path distance heatmap |",
        f"| `{OUT_EDGE_BAR.relative_to(ROOT)}` |  shortest-path edges  |",
        "",
        "## ",
        "",
        markdown_table(
            coverage_rows,
            ["system", "pkg", "ok", "missing_node", "no_path", "graph_error", "total_pairs_x_repeats", "ok_fraction"],
        ),
        "",
        "GetContacts  contact network key residues  `no_path`correlationplus/pytraj  communication path",
        "",
        "## BM213 vs C5apep focus pairs",
        "",
        "`delta_BM213_minus_C5apep < 0`  BM213  NetworkX shortest-path distance  pair  BM213  aggregate connection ",
        "",
        "### BM213  focus paths",
        "",
        markdown_table(
            bm_shorter,
            ["pkg", "pair_label", "apo", "BM213", "C5apep", "delta_BM213_minus_C5apep"],
            max_rows=12,
        ),
        "",
        "### C5apep  focus paths",
        "",
        markdown_table(
            c5_shorter,
            ["pkg", "pair_label", "apo", "BM213", "C5apep", "delta_BM213_minus_C5apep"],
            max_rows=12,
        ),
        "",
        "## Repeat-supported exact routes",
        "",
        " exact path  2  repeat  focus routes",
        "",
        markdown_table(
            repeated_focus,
            ["system", "pkg", "pair_label", "support_repeats", "replicates", "mean_path_distance", "path_residues"],
            max_rows=20,
        ),
        "",
        "## Consensus routes",
        "",
        " system/pkg  focus-pair route`consensus_intermediate_nodes_2of3`  2  repeat ",
        "",
        markdown_table(
            consensus_focus,
            [
                "system",
                "pkg",
                "pair_label",
                "successful_repeats",
                "mean_path_distance",
                "modal_path_support_repeats",
                "modal_path_residues",
                "consensus_intermediate_nodes_2of3",
                "consensus_edges_2of3",
            ],
            max_rows=30,
        ),
        "",
        "## ",
        "",
        markdown_table(
            top_nodes,
            ["system", "pkg", "node_residue", "node_context", "intermediate_occurrences", "unique_pairs", "focus_pair_occurrences"],
            max_rows=45,
        ),
        "",
        "##  shortest-path edges",
        "",
        markdown_table(
            top_edges,
            ["system", "pkg", "edge_label", "path_occurrences", "unique_pairs", "focus_pair_occurrences"],
            max_rows=45,
        ),
        "",
        "## ",
        "",
        "-  `shortest path proves activation pathway`",
        "-  `NetworkX shortest-path analysis of AlloViz-derived residue networks identified repeat-supported candidate communication routes connecting the ligand-proximal pocket, TM3 M120/R134, TM5 Y222, TM6 CWxP residues W255/Y258, and TM7/NPxxY-adjacent residues.`",
        "-  BM213 vs C5apep `networkx_focus_pair_distance_comparison.csv`  heatmapBM213  pair  stronger/G protein-biased  TM3-TM6/TM7C5apep  pair  partial efficacy ",
        "-  heatmap  PyMOL/PNG ",
        "",
        "## ",
        "",
        f"- Graph status rows: {len(status)}",
        f"- Key residues: {', '.join(KEY_RESIDUES_ORDER)}",
        f"- Focus pairs: {', '.join([f'{a}-{b}' for a, b in FOCUS_PAIRS])}",
        "",
    ]
    OUT_REPORT.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    TABLES_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    LOGS_DIR.mkdir(parents=True, exist_ok=True)

    raw_to_canonical, canonical_info, _key_to_canonical = load_mapping()
    paths, graph_status = compute_shortest_paths()
    summaries = summarize(paths, canonical_info)

    paths.to_csv(OUT_PATHS, index=False)
    summaries["exact_path_freq"].to_csv(OUT_PATH_FREQ, index=False)
    summaries["node_freq"].to_csv(OUT_NODE_FREQ, index=False)
    summaries["edge_freq"].to_csv(OUT_EDGE_FREQ, index=False)
    summaries["global_node_freq"].to_csv(OUT_GLOBAL_NODE_FREQ, index=False)
    summaries["global_edge_freq"].to_csv(OUT_GLOBAL_EDGE_FREQ, index=False)
    summaries["consensus"].to_csv(OUT_CONSENSUS, index=False)
    summaries["focus_compare"].to_csv(OUT_FOCUS_COMPARE, index=False)
    pd.DataFrame(graph_status).to_json(OUT_STATUS, orient="records", indent=2)

    plot_focus_heatmap(summaries["focus_compare"])
    plot_edge_support(summaries["global_edge_freq"])
    write_report(paths, summaries, graph_status)

    print(f"wrote {OUT_PATHS}")
    print(f"wrote {OUT_CONSENSUS}")
    print(f"wrote {OUT_REPORT}")
    print(f"wrote {OUT_HEATMAP}")


if __name__ == "__main__":
    main()
