#!/usr/bin/env python3
from __future__ import annotations

import html
import json
import re
from pathlib import Path

import pandas as pd


ROOT = Path("/mnt/e/work/modeling/c5ar/md")
OUTROOT = ROOT / "article/analysis/alloviz"
TABLES = OUTROOT / "tables"
HTML_DIR = OUTROOT / "figures/standalone_html"

SYSTEM_PDB = {
    "apo": ROOT / "article/analysis/alloviz/figures/pymol/apo_md-initial.pdb",
    "BM213": ROOT / "article/analysis/alloviz/figures/pymol/BM213_md-initial.pdb",
    "C5apep": ROOT / "article/analysis/alloviz/figures/pymol/C5apep_md-initial.pdb",
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


def js_string(text: str) -> str:
    return json.dumps(text)


def css_color(delta: float) -> str:
    if delta > 0:
        return "#1bb3b8"
    if delta < 0:
        return "#f28c1f"
    return "#e6e6e6"


def parse_ca_coords(pdb_path: Path) -> dict[int, dict[str, float]]:
    coords: dict[int, dict[str, float]] = {}
    for line in pdb_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        atom = line[12:16].strip()
        if atom != "CA":
            continue
        try:
            resid = int(line[22:26])
            coords[resid] = {
                "x": float(line[30:38]),
                "y": float(line[38:46]),
                "z": float(line[46:54]),
            }
        except ValueError:
            continue
    return coords


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


def make_nodes(key_nodes: pd.DataFrame, raw_map: dict[int, int], delta_col: str, coords: dict[int, dict[str, float]]) -> list[dict[str, object]]:
    nodes: list[dict[str, object]] = []
    sub = key_nodes[(key_nodes["pkg"] == "GetContacts") & key_nodes["biological_residue"].notna()]
    for row in sub.itertuples(index=False):
        delta = getattr(row, delta_col)
        if pd.isna(delta) or pd.isna(row.canonical_resid):
            continue
        raw = raw_map.get(int(row.canonical_resid))
        if raw is None or raw not in coords:
            continue
        delta = float(delta)
        nodes.append(
            {
                "label": str(row.biological_residue),
                "raw_resid": raw,
                "delta": delta,
                "color": css_color(delta),
                "radius": min(1.3, 0.55 + abs(delta) * 1.6),
                "coord": coords[raw],
            }
        )
    return nodes


def make_edges(
    edges: pd.DataFrame,
    raw_map: dict[int, int],
    labels: dict[int, str],
    delta_col: str,
    coords: dict[int, dict[str, float]],
    top_n: int = 15,
) -> list[dict[str, object]]:
    sub = edges[edges["pkg"] == "GetContacts"].copy()
    sub = sub.dropna(subset=[delta_col, "canonical_resid_a", "canonical_resid_b"])
    sub["abs_delta"] = sub[delta_col].abs()
    sub = sub.sort_values("abs_delta", ascending=False).head(top_n)
    out: list[dict[str, object]] = []
    for row in sub.itertuples(index=False):
        a = int(row.canonical_resid_a)
        b = int(row.canonical_resid_b)
        raw_a = raw_map.get(a)
        raw_b = raw_map.get(b)
        if raw_a is None or raw_b is None or raw_a not in coords or raw_b not in coords:
            continue
        delta = float(getattr(row, delta_col))
        out.append(
            {
                "label": f"{labels.get(a, a)}-{labels.get(b, b)}",
                "delta": delta,
                "color": css_color(delta),
                "radius": min(0.45, 0.12 + abs(delta) * 0.65),
                "start": coords[raw_a],
                "end": coords[raw_b],
            }
        )
    return out


def html_doc(title: str, pdb_text: str, nodes: list[dict[str, object]], edges: list[dict[str, object]], positive: str, negative: str) -> str:
    return f"""<!doctype html>
<html>
<head>
  <meta charset="utf-8">
  <title>{html.escape(title)}</title>
  <script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
  <style>
    body {{ margin: 0; font-family: Arial, sans-serif; }}
    #viewer {{ width: 100vw; height: 100vh; }}
    #legend {{
      position: absolute; left: 14px; top: 12px; z-index: 10;
      background: rgba(255,255,255,0.90); padding: 10px 12px;
      border: 1px solid #ccc; border-radius: 4px; font-size: 13px;
      max-width: 360px;
    }}
    .swatch {{ display:inline-block; width:12px; height:12px; margin-right:6px; vertical-align:middle; }}
  </style>
</head>
<body>
  <div id="legend">
    <b>{html.escape(title)}</b><br>
    <span class="swatch" style="background:#1bb3b8"></span>{html.escape(positive)}<br>
    <span class="swatch" style="background:#f28c1f"></span>{html.escape(negative)}<br>
    spheres/sticks = key residues, cylinders = top delta edges
  </div>
  <div id="viewer"></div>
  <script>
    const pdb = {js_string(pdb_text)};
    const nodes = {json.dumps(nodes)};
    const edges = {json.dumps(edges)};
    const viewer = $3Dmol.createViewer("viewer", {{backgroundColor: "white"}});
    viewer.addModel(pdb, "pdb");
    viewer.setStyle({{}}, {{cartoon: {{color: "lightgray", opacity: 0.72}}}});
    for (const node of nodes) {{
      viewer.setStyle({{resi: node.raw_resid}}, {{
        stick: {{colorscheme: "grayCarbon", radius: 0.17}},
        sphere: {{radius: node.radius, color: node.color}}
      }});
      viewer.addLabel(node.label, {{
        position: node.coord,
        fontSize: 13,
        fontColor: "black",
        backgroundColor: "white",
        backgroundOpacity: 0.65,
        inFront: true
      }});
    }}
    for (const edge of edges) {{
      viewer.addCylinder({{
        start: edge.start,
        end: edge.end,
        radius: edge.radius,
        color: edge.color,
        opacity: 0.85,
        dashed: true
      }});
    }}
    viewer.zoomTo();
    viewer.render();
  </script>
</body>
</html>
"""


def write_readme(rows: list[dict[str, str]]) -> None:
    lines = [
        "# Standalone browser HTML views",
        "",
        " HTML Jupyter notebook / nglview widget ",
        "",
        " 3Dmol.js CDN `https://3Dmol.csb.pitt.edu` `../pymol/*.pml`",
        "",
        "| comparison | file | display structure |",
        "|---|---|---|",
    ]
    for row in rows:
        lines.append(f"| {row['comparison']} | `{row['file']}` | {row['display_system']} |")
    (HTML_DIR / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    HTML_DIR.mkdir(parents=True, exist_ok=True)
    mapping = pd.read_csv(TABLES / "residue_mapping_alloviz.csv")
    key_nodes = pd.read_csv(TABLES / "key_residue_network_summary.csv")
    edge_tables = {
        "BM213_vs_apo": pd.read_csv(TABLES / "delta_edges_BM213_vs_apo.csv"),
        "C5apep_vs_apo": pd.read_csv(TABLES / "delta_edges_C5apep_vs_apo.csv"),
        "BM213_vs_C5apep": pd.read_csv(TABLES / "delta_edges_BM213_vs_C5apep.csv"),
    }

    rows: list[dict[str, str]] = []
    for comp in COMPARISONS:
        display = comp["display_system"]
        pdb_path = SYSTEM_PDB[display]
        pdb_text = pdb_path.read_text(encoding="utf-8", errors="replace")
        coords = parse_ca_coords(pdb_path)
        raw_map = raw_resid_map(mapping, display)
        labels = label_map(mapping, display)
        nodes = make_nodes(key_nodes, raw_map, comp["delta_col"], coords)
        edges = make_edges(edge_tables[comp["name"]], raw_map, labels, comp["delta_col"], coords)
        title = f"{comp['name']} GetContacts network"
        out = HTML_DIR / f"{comp['name']}_GetContacts_network_standalone.html"
        out.write_text(
            html_doc(title, pdb_text, nodes, edges, comp["positive_label"], comp["negative_label"]),
            encoding="utf-8",
        )
        rows.append(
            {
                "comparison": comp["name"],
                "file": out.name,
                "display_system": display,
            }
        )
    write_readme(rows)
    print(f"Wrote standalone HTML under {HTML_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

