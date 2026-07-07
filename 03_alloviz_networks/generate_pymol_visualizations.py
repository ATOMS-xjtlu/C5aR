#!/usr/bin/env python3
from __future__ import annotations

import re
from pathlib import Path

import pandas as pd


ROOT = Path("/mnt/e/work/modeling/c5ar/md")
OUTROOT = ROOT / "article/analysis/alloviz"
TABLES = OUTROOT / "tables"
PYMOL_DIR = OUTROOT / "figures/pymol"
PYMOL_DIR_WINDOWS = "E:/work/modeling/c5ar/md/article/analysis/alloviz/figures/pymol"

SYSTEM_PDB = {
    "apo": "article/analysis/alloviz/figures/pymol/apo_md-initial.pdb",
    "BM213": "article/analysis/alloviz/figures/pymol/BM213_md-initial.pdb",
    "C5apep": "article/analysis/alloviz/figures/pymol/C5apep_md-initial.pdb",
}

COMPARISONS = [
    {
        "name": "BM213_vs_apo",
        "display_system": "BM213",
        "delta_col": "delta_BM213_minus_apo",
        "positive_label": "BM213 higher",
        "negative_label": "apo higher",
    },
    {
        "name": "C5apep_vs_apo",
        "display_system": "C5apep",
        "delta_col": "delta_C5apep_minus_apo",
        "positive_label": "C5apep higher",
        "negative_label": "apo higher",
    },
    {
        "name": "BM213_vs_C5apep",
        "display_system": "BM213",
        "delta_col": "delta_BM213_minus_C5apep",
        "positive_label": "BM213 higher",
        "negative_label": "C5apep higher",
    },
]


def safe_name(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9_]+", "_", str(text)).strip("_")


def load_mapping() -> pd.DataFrame:
    mapping = pd.read_csv(TABLES / "residue_mapping_alloviz.csv")
    return mapping


def raw_resid_map(mapping: pd.DataFrame, system: str) -> dict[int, int]:
    sub = mapping[mapping["system"] == system]
    return {
        int(row.canonical_resid): int(row.raw_resid)
        for row in sub.itertuples(index=False)
        if pd.notna(row.canonical_resid) and pd.notna(row.raw_resid)
    }


def label_map(mapping: pd.DataFrame, system: str) -> dict[int, str]:
    sub = mapping[mapping["system"] == system]
    labels: dict[int, str] = {}
    for row in sub.itertuples(index=False):
        if pd.isna(row.canonical_resid):
            continue
        label = row.biological_residue if isinstance(row.biological_residue, str) else row.canonical_label
        labels[int(row.canonical_resid)] = str(label)
    return labels


def pml_header(display_system: str, title: str) -> list[str]:
    pdb = f"{PYMOL_DIR_WINDOWS}/{Path(SYSTEM_PDB[display_system]).name}"
    return [
        f"# {title}",
        "# Open in Windows PyMOL with File > Run Script",
        "reinitialize",
        f'load "{pdb}", receptor',
        "remove receptor and not polymer.protein",
        "hide everything",
        "bg_color white",
        "set ray_opaque_background, off",
        "set antialias, 2",
        "set ambient, 0.35",
        "set cartoon_fancy_helices, on",
        "show cartoon, receptor",
        "color gray85, receptor",
        "set cartoon_transparency, 0.20, receptor",
        "set_color delta_positive, [0.10, 0.70, 0.72]",
        "set_color delta_negative, [0.95, 0.55, 0.12]",
        "set_color delta_neutral, [0.95, 0.95, 0.95]",
        "set_color edge_positive, [0.00, 0.55, 0.60]",
        "set_color edge_negative, [0.90, 0.42, 0.00]",
        "set_color key_gray, [0.25, 0.25, 0.25]",
        "set label_size, 16",
        "set label_color, black",
        "set dash_gap, 0.20",
        "set dash_radius, 0.07",
        "set dash_round_ends, on",
        "",
    ]


def add_node_commands(
    lines: list[str],
    nodes: pd.DataFrame,
    display_system: str,
    raw_map: dict[int, int],
    delta_col: str,
) -> None:
    lines.extend(
        [
            "# Key residue nodes. Cyan means positive state higher. Orange means comparator higher.",
            "select all_key_nodes, none",
        ]
    )
    for row in nodes.itertuples(index=False):
        if pd.isna(getattr(row, delta_col)) or pd.isna(row.canonical_resid):
            continue
        canonical = int(row.canonical_resid)
        raw = raw_map.get(canonical)
        if raw is None:
            continue
        label = row.biological_residue if isinstance(row.biological_residue, str) else row.canonical_label
        delta = float(getattr(row, delta_col))
        color = "delta_positive" if delta > 0 else "delta_negative" if delta < 0 else "delta_neutral"
        scale = min(0.75, 0.25 + abs(delta) * 1.4)
        sel = safe_name(f"node_{label}_{raw}")
        lines.extend(
            [
                f"select {sel}, receptor and resi {raw}",
                f"show sticks, {sel}",
                f"color {color}, {sel}",
                f"show spheres, {sel} and name CA",
                f"set sphere_scale, {scale:.3f}, {sel} and name CA",
                f'label {sel} and name CA, "{label}"',
                f"select all_key_nodes, all_key_nodes or {sel}",
            ]
        )
    lines.extend(["", "orient all_key_nodes", "zoom all_key_nodes, 12", ""])


def add_edge_commands(
    lines: list[str],
    edges: pd.DataFrame,
    raw_map: dict[int, int],
    labels: dict[int, str],
    delta_col: str,
    top_n: int = 15,
) -> None:
    edges = edges[edges["pkg"] == "GetContacts"].copy()
    edges = edges.dropna(subset=[delta_col, "canonical_resid_a", "canonical_resid_b"])
    edges["abs_delta"] = edges[delta_col].abs()
    edges = edges.sort_values("abs_delta", ascending=False).head(top_n)
    lines.extend(["# Top delta edges shown as CA-CA dashed links.", "select all_delta_edges, none"])
    for idx, row in enumerate(edges.itertuples(index=False), start=1):
        a = int(row.canonical_resid_a)
        b = int(row.canonical_resid_b)
        raw_a = raw_map.get(a)
        raw_b = raw_map.get(b)
        if raw_a is None or raw_b is None:
            continue
        delta = float(getattr(row, delta_col))
        color = "edge_positive" if delta > 0 else "edge_negative"
        label_a = labels.get(a, f"{a}")
        label_b = labels.get(b, f"{b}")
        edge_name = safe_name(f"edge_{idx}_{label_a}_{label_b}")
        width = min(7.0, 2.0 + abs(delta) * 8.0)
        lines.extend(
            [
                f"distance {edge_name}, receptor and resi {raw_a} and name CA, receptor and resi {raw_b} and name CA",
                f"hide labels, {edge_name}",
                f"color {color}, {edge_name}",
                f"set dash_width, {width:.2f}, {edge_name}",
                f"select all_delta_edges, all_delta_edges or receptor and resi {raw_a}+{raw_b}",
            ]
        )
    lines.extend(["", ""])


def write_index(rows: list[dict[str, str]]) -> None:
    lines = [
        "# PyMOL visualization outputs",
        "",
        " PyMOL-ready  ray-traced  PNG/PSE WSL  `pymol` ",
        "",
        " Windows PyMOL  `.pml` `pymol <script.pml>` `ray 2400,1800`  `png <name>.png, dpi=300`",
        "",
        "| comparison | display structure | script | meaning |",
        "|---|---|---|---|",
    ]
    for row in rows:
        lines.append(
            f"| {row['comparison']} | {row['display_system']} | `{row['script']}` | {row['meaning']} |"
        )
    lines.extend(
        [
            "",
            "",
            "",
            "- cyan/turquoise `BM213_vs_C5apep`  BM213 ",
            "- orange `BM213_vs_C5apep`  C5apep ",
            "- residue sticks/spheres",
            "- dashed linesGetContacts top delta edges CA-CA",
            "",
            " `tables/`  repeat-aware CSV ",
        ]
    )
    (PYMOL_DIR / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    PYMOL_DIR.mkdir(parents=True, exist_ok=True)
    mapping = load_mapping()
    key_nodes = pd.read_csv(TABLES / "key_residue_network_summary.csv")
    edge_tables = {
        "BM213_vs_apo": pd.read_csv(TABLES / "delta_edges_BM213_vs_apo.csv"),
        "C5apep_vs_apo": pd.read_csv(TABLES / "delta_edges_C5apep_vs_apo.csv"),
        "BM213_vs_C5apep": pd.read_csv(TABLES / "delta_edges_BM213_vs_C5apep.csv"),
    }

    rows: list[dict[str, str]] = []
    for comp in COMPARISONS:
        name = comp["name"]
        display = comp["display_system"]
        delta_col = comp["delta_col"]
        raw_map = raw_resid_map(mapping, display)
        labels = label_map(mapping, display)

        nodes = key_nodes[
            (key_nodes["pkg"] == "GetContacts") & key_nodes["biological_residue"].notna()
        ].copy()
        lines = pml_header(display, f"{name} GetContacts node/edge delta")
        lines.extend(
            [
                f"# Positive color: {comp['positive_label']}",
                f"# Negative color: {comp['negative_label']}",
                "",
            ]
        )
        add_node_commands(lines, nodes, display, raw_map, delta_col)
        add_edge_commands(lines, edge_tables[name], raw_map, labels, delta_col, top_n=15)
        lines.extend(
            [
                "set stick_radius, 0.16, all_key_nodes",
                "show cartoon, receptor",
                "zoom all_key_nodes, 10",
                "# Optional rendering commands:",
                "# ray 2400, 1800",
                f"# png {name}_GetContacts_network.png, dpi=300",
                "",
            ]
        )

        script = PYMOL_DIR / f"{name}_GetContacts_network.pml"
        script.write_text("\n".join(lines) + "\n", encoding="utf-8")
        rows.append(
            {
                "comparison": name,
                "display_system": display,
                "script": script.relative_to(OUTROOT).as_posix(),
                "meaning": f"{comp['positive_label']} vs {comp['negative_label']}",
            }
        )

    write_index(rows)
    print(f"Wrote PyMOL scripts under {PYMOL_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
