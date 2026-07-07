#!/usr/bin/env python3
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from generate_standalone_3dmol_html import (  # noqa: E402
    COMPARISONS,
    OUTROOT,
    SYSTEM_PDB,
    TABLES,
    label_map,
    make_edges,
    make_nodes,
    parse_ca_coords,
    raw_resid_map,
)


PNG_DIR = OUTROOT / "figures/png"


def coord_array(coords: dict[int, dict[str, float]], residues: list[int]) -> np.ndarray:
    pts = []
    for resid in residues:
        if resid in coords:
            c = coords[resid]
            pts.append([c["x"], c["y"], c["z"]])
    return np.asarray(pts, dtype=float)


def set_axes_equal(ax, pts: np.ndarray) -> None:
    mins = pts.min(axis=0)
    maxs = pts.max(axis=0)
    center = (mins + maxs) / 2
    radius = (maxs - mins).max() / 2
    ax.set_xlim(center[0] - radius, center[0] + radius)
    ax.set_ylim(center[1] - radius, center[1] + radius)
    ax.set_zlim(center[2] - radius, center[2] + radius)


def plot_network(comp: dict[str, str], mapping: pd.DataFrame, key_nodes: pd.DataFrame, edge_tables: dict[str, pd.DataFrame]) -> Path:
    display = comp["display_system"]
    pdb_path = SYSTEM_PDB[display]
    coords = parse_ca_coords(pdb_path)
    raw_map = raw_resid_map(mapping, display)
    labels = label_map(mapping, display)
    nodes = make_nodes(key_nodes, raw_map, comp["delta_col"], coords)
    edges = make_edges(edge_tables[comp["name"]], raw_map, labels, comp["delta_col"], coords, top_n=15)

    receptor_resids = sorted(set(raw_map.values()))
    backbone = coord_array(coords, receptor_resids)
    if backbone.size == 0:
        raise ValueError(f"No CA coordinates found for {comp['name']}")

    fig = plt.figure(figsize=(8, 8), dpi=300)
    ax = fig.add_subplot(111, projection="3d")
    ax.plot(backbone[:, 0], backbone[:, 1], backbone[:, 2], color="#d0d0d0", lw=1.0, alpha=0.78)

    for edge in edges:
        start = edge["start"]
        end = edge["end"]
        delta = float(edge["delta"])
        lw = min(7.0, 1.8 + abs(delta) * 9.0)
        ax.plot(
            [start["x"], end["x"]],
            [start["y"], end["y"]],
            [start["z"], end["z"]],
            color=edge["color"],
            lw=lw,
            alpha=0.78,
        )

    for node in nodes:
        coord = node["coord"]
        delta = float(node["delta"])
        size = min(520, 110 + abs(delta) * 1600)
        ax.scatter(
            [coord["x"]],
            [coord["y"]],
            [coord["z"]],
            s=size,
            color=node["color"],
            edgecolor="black",
            linewidth=0.45,
            depthshade=False,
            zorder=8,
        )
        ax.text(
            coord["x"] + 0.9,
            coord["y"] + 0.9,
            coord["z"] + 0.9,
            str(node["label"]),
            fontsize=8.0,
            color="black",
            zorder=9,
        )

    title = f"{comp['name']} GetContacts network"
    ax.set_title(title, fontsize=13, pad=12)
    ax.text2D(
        0.02,
        0.96,
        f"cyan: {comp['positive_label']}\norange: {comp['negative_label']}\nNodes: key residues | Lines: top delta edges",
        transform=ax.transAxes,
        fontsize=8,
        va="top",
        bbox=dict(facecolor="white", edgecolor="#cccccc", alpha=0.88, boxstyle="round,pad=0.35"),
    )

    set_axes_equal(ax, backbone)
    ax.view_init(elev=15, azim=-72)
    ax.set_axis_off()
    fig.tight_layout()

    out = PNG_DIR / f"{comp['name']}_GetContacts_network.png"
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out


def write_readme(outputs: list[Path]) -> None:
    lines = [
        "# Static PNG GetContacts Network Views",
        "",
        " PNG  HTML  GetContacts node/edge  PDB CA ",
        "",
        "| file |",
        "|---|",
    ]
    for out in outputs:
        lines.append(f"| `{out.name}` |")
    lines.extend(
        [
            "",
            "",
            "",
            "- cyan/turquoise",
            "- orange",
            "- receptor CA backbone trace",
            "- ",
            "- GetContacts top delta edges",
        ]
    )
    (PNG_DIR / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    PNG_DIR.mkdir(parents=True, exist_ok=True)
    mapping = pd.read_csv(TABLES / "residue_mapping_alloviz.csv")
    key_nodes = pd.read_csv(TABLES / "key_residue_network_summary.csv")
    edge_tables = {
        "BM213_vs_apo": pd.read_csv(TABLES / "delta_edges_BM213_vs_apo.csv"),
        "C5apep_vs_apo": pd.read_csv(TABLES / "delta_edges_C5apep_vs_apo.csv"),
        "BM213_vs_C5apep": pd.read_csv(TABLES / "delta_edges_BM213_vs_C5apep.csv"),
    }
    outputs = [plot_network(comp, mapping, key_nodes, edge_tables) for comp in COMPARISONS]
    write_readme(outputs)
    for out in outputs:
        print(out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
