#!/usr/bin/env python3
from __future__ import annotations

import html
import json
import os
import re
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-c5ar-alloviz")

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import pandas as pd


ROOT = Path("/mnt/e/work/modeling/c5ar/md")
OUTROOT = ROOT / "article/analysis/alloviz"
TABLES = OUTROOT / "tables"
ACTIVE_DIR = OUTROOT / "official_views/nglview/active"
EDGE_SCORES = TABLES / "edge_scores_long.csv"

SYSTEM_PDB = {
    "BM213": OUTROOT / "figures/pymol/BM213_md-initial.pdb",
    "C5apep": OUTROOT / "figures/pymol/C5apep_md-initial.pdb",
}

SYSTEM_PDB_WINDOWS = {
    "BM213": "E:/work/modeling/c5ar/md/article/analysis/alloviz/figures/pymol/BM213_md-initial.pdb",
    "C5apep": "E:/work/modeling/c5ar/md/article/analysis/alloviz/figures/pymol/C5apep_md-initial.pdb",
}


@dataclass(frozen=True)
class RouteSpec:
    system: str
    pkg: str
    pair_label: str
    route_label: str
    class_name: str


SYSTEM_STYLES = {
    "BM213": {
        "title": "BM213 candidate activation-like communication routes",
        "slug": "BM213_activation_route",
        "main_color": "#0fa7ad",
        "intermediate_color": "#9bdfe1",
        "supplement_color": "#4866c8",
        "pymol_main": "[0.05, 0.65, 0.68]",
        "pymol_intermediate": "[0.55, 0.86, 0.88]",
        "pymol_supplement": "[0.28, 0.40, 0.78]",
    },
    "C5apep": {
        "title": "C5apep candidate activation-like communication routes",
        "slug": "C5apep_activation_route",
        "main_color": "#f28c1f",
        "intermediate_color": "#ffd39a",
        "supplement_color": "#b35a00",
        "pymol_main": "[0.95, 0.55, 0.12]",
        "pymol_intermediate": "[1.00, 0.78, 0.48]",
        "pymol_supplement": "[0.70, 0.35, 0.00]",
    },
}

ROUTES = [
    RouteSpec("BM213", "GetContacts", "M120-W255", "TM3-CWxP direct bridge", "main"),
    RouteSpec("BM213", "GetContacts", "W255-Y258", "CWxP microswitch core", "main"),
    RouteSpec("BM213", "GetContacts", "W255-Y290", "TM6-to-TM7/IWI route", "main"),
    RouteSpec("BM213", "GetContacts", "R134-Y222", "DRY-like-to-TM5 route", "main"),
    RouteSpec("BM213", "GetContacts", "S171-Y258", "pocket-to-CWxP route", "main"),
    RouteSpec("BM213", "correlationplus_CA_Pear", "I116-Y258", "IWI-to-CWxP correlation route", "supplement"),
    RouteSpec("BM213", "correlationplus_CA_Pear", "M120-Y258", "TM3-to-CWxP correlation route", "supplement"),
    RouteSpec("BM213", "correlationplus_CA_Pear", "M120-W255", "TM3-CWxP correlation route", "supplement"),
    RouteSpec("C5apep", "GetContacts", "F251-W255", "PIF/CWxP local route", "main"),
    RouteSpec("C5apep", "GetContacts", "N296-Y300", "TM7/NPxxY-adjacent contact route", "main"),
    RouteSpec("C5apep", "GetContacts", "I91-M120", "IWI-to-TM3/CWxP contact route", "supplement"),
    RouteSpec("C5apep", "correlationplus_CA_Pear", "M120-Y300", "TM3-to-NPxxY correlation route", "supplement"),
    RouteSpec("C5apep", "correlationplus_CA_Pear", "R134-Y222", "DRY-like-to-TM5 correlation route", "supplement"),
]

EXTRA_ROUTES = {
    "BM213": [
        {
            "pkg": "GetContacts_bridge",
            "pair_label": "W255-L127-Y222",
            "route_label": "W255-L127-Y222 bridge",
            "class_name": "supplement",
            "path_labels": ["W255", "L127", "Y222"],
            "support_repeats": 2,
            "successful_repeats": 3,
            "mean_path_distance": None,
        },
        {
            "pkg": "GetContacts_bridge",
            "pair_label": "W255-L127-R134",
            "route_label": "W255-L127-R134 bridge",
            "class_name": "supplement",
            "path_labels": ["W255", "L127", "R134"],
            "support_repeats": 2,
            "successful_repeats": 3,
            "mean_path_distance": None,
        },
        {
            "pkg": "manual_link",
            "pair_label": "Y222-F251",
            "route_label": "Y222-F251 link",
            "class_name": "main",
            "path_labels": ["Y222", "F251"],
            "support_repeats": 1,
            "successful_repeats": 3,
            "mean_path_distance": None,
        },
    ],
    "C5apep": [
        {
            "pkg": "manual_link",
            "pair_label": "S171-R206",
            "route_label": "peptide-contact route",
            "class_name": "main",
            "path_labels": ["S171", "R206"],
            "support_repeats": 3,
            "successful_repeats": 3,
            "mean_path_distance": None,
        },
        {
            "pkg": "manual_link",
            "pair_label": "Y222-F251",
            "route_label": "Y222-F251 link",
            "class_name": "main",
            "path_labels": ["Y222", "F251"],
            "support_repeats": 1,
            "successful_repeats": 3,
            "mean_path_distance": None,
        },
    ]
}

EXTRA_NODES: dict[str, list[dict[str, object]]] = {}


def safe_name(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9_]+", "_", text).strip("_")


def js_string(text: str) -> str:
    return json.dumps(text)


def parse_ca_coords(pdb_path: Path) -> dict[int, dict[str, float]]:
    coords: dict[int, dict[str, float]] = {}
    for line in pdb_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        if line[12:16].strip() != "CA":
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


def label_to_canonical_map(mapping: pd.DataFrame, system: str) -> dict[str, int]:
    sub = mapping[mapping["system"] == system]
    out: dict[str, int] = {}
    priorities: dict[str, int] = {}

    def assign(label: object, canonical: int, priority: int) -> None:
        if not isinstance(label, str) or not label:
            return
        if priority >= priorities.get(label, -1):
            out[label] = canonical
            priorities[label] = priority

    for row in sub.itertuples(index=False):
        if pd.isna(row.canonical_resid):
            continue
        canonical = int(row.canonical_resid)
        has_biological_label = isinstance(row.biological_residue, str) and bool(row.biological_residue)
        is_key_residue = bool(getattr(row, "is_key_residue", False))
        assign(row.canonical_label, canonical, 2 if has_biological_label or is_key_residue else 0)
        assign(row.biological_residue, canonical, 3)
        if isinstance(row.aa, str):
            assign(f"{row.aa}{canonical + 33}", canonical, 1)
    return out


def canonical_label_map(mapping: pd.DataFrame, system: str) -> dict[int, str]:
    sub = mapping[mapping["system"] == system]
    labels: dict[int, str] = {}
    for row in sub.itertuples(index=False):
        if pd.isna(row.canonical_resid):
            continue
        canonical = int(row.canonical_resid)
        if isinstance(row.biological_residue, str):
            label = row.biological_residue
        elif isinstance(row.aa, str):
            label = f"{row.aa}{canonical + 33}"
        else:
            label = row.canonical_label
        labels[canonical] = str(label)
    return labels


def context_map(mapping: pd.DataFrame, system: str) -> dict[int, str]:
    sub = mapping[mapping["system"] == system]
    out: dict[int, str] = {}
    for row in sub.itertuples(index=False):
        if pd.isna(row.canonical_resid):
            continue
        context = row.panel_context if isinstance(row.panel_context, str) else ""
        out[int(row.canonical_resid)] = context
    return out


def gpcrdb_map(mapping: pd.DataFrame, system: str) -> dict[int, str]:
    sub = mapping[mapping["system"] == system]
    out: dict[int, str] = {}
    for row in sub.itertuples(index=False):
        if pd.isna(row.canonical_resid):
            continue
        gpcrdb = row.gpcrdb if isinstance(row.gpcrdb, str) else ""
        out[int(row.canonical_resid)] = gpcrdb
    return out


def parse_path(path_text: str) -> list[str]:
    return [part.strip() for part in str(path_text).split("->") if part.strip()]


def route_row(consensus: pd.DataFrame, spec: RouteSpec) -> pd.Series:
    sub = consensus[
        (consensus["system"] == spec.system)
        & (consensus["pkg"] == spec.pkg)
        & (consensus["pair_label"] == spec.pair_label)
    ]
    if sub.empty:
        raise ValueError(f"Missing route row: {spec.system} {spec.pkg} {spec.pair_label}")
    return sub.iloc[0]


def route_color(system: str, class_name: str) -> str:
    style = SYSTEM_STYLES[system]
    return style["supplement_color"] if class_name == "supplement" else style["main_color"]


def strength_pkg(pkg: str) -> str:
    if pkg in {"GetContacts_bridge", "manual_link"}:
        return "GetContacts"
    return pkg


def segment_strength(edge_scores: pd.DataFrame, system: str, pkg: str, a: int, b: int) -> tuple[float | None, int]:
    sub = edge_scores[
        (edge_scores["system"] == system)
        & (edge_scores["pkg"] == strength_pkg(pkg))
        & (
            (
                (edge_scores["canonical_resid_a"] == a)
                & (edge_scores["canonical_resid_b"] == b)
            )
            | (
                (edge_scores["canonical_resid_a"] == b)
                & (edge_scores["canonical_resid_b"] == a)
            )
        )
    ].copy()
    if sub.empty:
        return None, 0
    sub["abs_score"] = sub["score"].abs()
    positive = sub[sub["abs_score"] > 0]
    if positive.empty:
        return 0.0, 0
    return float(positive["abs_score"].mean()), int(positive["replicate"].nunique())


def annotate_route_segments(routes: list[dict[str, object]], edge_scores: pd.DataFrame) -> None:
    all_scores: dict[tuple[str, str], list[float]] = {}
    for route in routes:
        system = str(route["system"])
        pkg = str(route["pkg"])
        strength_method = strength_pkg(pkg)
        canonical_path = [int(item) for item in route["path_canonical_resids"]]
        labels = [str(item) for item in route["path_labels"]]
        raw_path = [int(item) for item in route["path_raw_resids"]]
        segments = []
        for index, (canonical_a, canonical_b, raw_a, raw_b, label_a, label_b) in enumerate(
            zip(canonical_path, canonical_path[1:], raw_path, raw_path[1:], labels, labels[1:]),
            start=1,
        ):
            score, support = segment_strength(edge_scores, system, pkg, canonical_a, canonical_b)
            if score is not None and score > 0:
                all_scores.setdefault((system, strength_method), []).append(score)
            segments.append(
                {
                    "index": index,
                    "canonical_a": canonical_a,
                    "canonical_b": canonical_b,
                    "raw_a": raw_a,
                    "raw_b": raw_b,
                    "label": f"{label_a}-{label_b}",
                    "strength_pkg": strength_method,
                    "strength_score": score,
                    "strength_support_repeats": support,
                }
            )
        route["segments"] = segments

    ranges: dict[tuple[str, str], tuple[float, float]] = {}
    for key, values in all_scores.items():
        ranges[key] = (min(values), max(values))

    for route in routes:
        system = str(route["system"])
        strength_method = strength_pkg(str(route["pkg"]))
        min_score, max_score = ranges.get((system, strength_method), (0.0, 1.0))
        for segment in route["segments"]:
            score = segment["strength_score"]
            if score is None or score <= 0:
                norm = 0.0
            elif max_score <= min_score:
                norm = 1.0
            else:
                norm = (float(score) - min_score) / (max_score - min_score)
                norm = max(0.0, min(1.0, norm))
            segment["strength_norm"] = norm
            segment["pymol_dash_width"] = round(3.0 + 6.0 * norm, 2)
            segment["pymol_dash_radius"] = round(0.06 + 0.12 * norm, 3)
            segment["html_radius"] = round(0.14 + 0.30 * norm, 3)


def build_routes(
    system: str,
    consensus: pd.DataFrame,
    mapping: pd.DataFrame,
    coords: dict[int, dict[str, float]],
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[str]]:
    raw_map = raw_resid_map(mapping, system)
    label_to_canonical = label_to_canonical_map(mapping, system)
    labels = canonical_label_map(mapping, system)
    contexts = context_map(mapping, system)
    gpcrdb = gpcrdb_map(mapping, system)
    warnings: list[str] = []
    routes: list[dict[str, object]] = []
    node_index: dict[int, dict[str, object]] = {}

    for spec in [item for item in ROUTES if item.system == system]:
        row = route_row(consensus, spec)
        labels_in_path = parse_path(row.modal_path_residues)
        canonical_path: list[int] = []
        raw_path: list[int] = []
        missing: list[str] = []
        for label in labels_in_path:
            canonical = label_to_canonical.get(label)
            if canonical is None:
                missing.append(label)
                continue
            raw = raw_map.get(canonical)
            if raw is None or raw not in coords:
                missing.append(label)
                continue
            canonical_path.append(canonical)
            raw_path.append(raw)
            is_endpoint = label in {labels_in_path[0], labels_in_path[-1]}
            old = node_index.get(canonical)
            if old is None:
                node_index[canonical] = {
                    "canonical_resid": canonical,
                    "raw_resid": raw,
                    "label": labels.get(canonical, label),
                    "context": contexts.get(canonical, ""),
                    "gpcrdb": gpcrdb.get(canonical, ""),
                    "coord": coords[raw],
                    "endpoint": is_endpoint,
                    "classes": {spec.class_name},
                }
            else:
                old["endpoint"] = bool(old["endpoint"] or is_endpoint)
                old["classes"].add(spec.class_name)
        if missing:
            warnings.append(f"{system} {spec.pkg} {spec.pair_label}: missing {', '.join(missing)}")
        if len(raw_path) < 2:
            warnings.append(f"{system} {spec.pkg} {spec.pair_label}: fewer than two mappable residues")
            continue
        support = int(row.modal_path_support_repeats) if pd.notna(row.modal_path_support_repeats) else 0
        distance = float(row.mean_path_distance) if pd.notna(row.mean_path_distance) else None
        routes.append(
            {
                "system": system,
                "pkg": spec.pkg,
                "pair_label": spec.pair_label,
                "route_label": spec.route_label,
                "class_name": spec.class_name,
                "path_labels": [labels.get(c, str(c)) for c in canonical_path],
                "path_canonical_resids": canonical_path,
                "path_raw_resids": raw_path,
                "support_repeats": support,
                "successful_repeats": int(row.successful_repeats),
                "mean_path_distance": distance,
                "color": route_color(system, spec.class_name),
                "radius": 0.26 if spec.class_name == "main" else 0.14,
                "opacity": 0.82 if spec.class_name == "main" else 0.55,
            }
        )

    for extra in EXTRA_ROUTES.get(system, []):
        canonical_path: list[int] = []
        raw_path: list[int] = []
        missing: list[str] = []
        extra_labels = list(extra["path_labels"])
        for label in extra_labels:
            canonical = label_to_canonical.get(label)
            if canonical is None:
                missing.append(label)
                continue
            raw = raw_map.get(canonical)
            if raw is None or raw not in coords:
                missing.append(label)
                continue
            canonical_path.append(canonical)
            raw_path.append(raw)
            old = node_index.get(canonical)
            is_endpoint = label in {extra_labels[0], extra_labels[-1]}
            if old is None:
                node_index[canonical] = {
                    "canonical_resid": canonical,
                    "raw_resid": raw,
                    "label": labels.get(canonical, label),
                    "context": contexts.get(canonical, ""),
                    "gpcrdb": gpcrdb.get(canonical, ""),
                    "coord": coords[raw],
                    "endpoint": is_endpoint,
                    "classes": {str(extra["class_name"])},
                }
            else:
                old["endpoint"] = bool(old["endpoint"] or is_endpoint)
                old["classes"].add(str(extra["class_name"]))
        if missing:
            warnings.append(f"{system} {extra['pkg']} {extra['pair_label']}: missing {', '.join(missing)}")
        if len(raw_path) < 2:
            warnings.append(f"{system} {extra['pkg']} {extra['pair_label']}: fewer than two mappable residues")
            continue
        routes.append(
            {
                "system": system,
                "pkg": extra["pkg"],
                "pair_label": extra["pair_label"],
                "route_label": extra["route_label"],
                "class_name": extra["class_name"],
                "path_labels": [labels.get(c, str(c)) for c in canonical_path],
                "path_canonical_resids": canonical_path,
                "path_raw_resids": raw_path,
                "support_repeats": int(extra["support_repeats"]),
                "successful_repeats": int(extra["successful_repeats"]),
                "mean_path_distance": extra["mean_path_distance"],
                "color": route_color(system, str(extra["class_name"])),
                "radius": 0.14,
                "opacity": 0.62,
            }
        )

    for extra_node in EXTRA_NODES.get(system, []):
        label = str(extra_node["label"])
        canonical = label_to_canonical.get(label)
        if canonical is None:
            warnings.append(f"{system} extra node {label}: missing")
            continue
        raw = raw_map.get(canonical)
        if raw is None or raw not in coords:
            warnings.append(f"{system} extra node {label}: missing coordinates")
            continue
        old = node_index.get(canonical)
        class_name = str(extra_node["class_name"])
        is_endpoint = bool(extra_node.get("endpoint", True))
        if old is None:
            node_index[canonical] = {
                "canonical_resid": canonical,
                "raw_resid": raw,
                "label": labels.get(canonical, label),
                "context": contexts.get(canonical, ""),
                "gpcrdb": gpcrdb.get(canonical, ""),
                "coord": coords[raw],
                "endpoint": is_endpoint,
                "classes": {class_name},
            }
        else:
            old["endpoint"] = bool(old["endpoint"] or is_endpoint)
            old["classes"].add(class_name)

    nodes = []
    for node in node_index.values():
        classes = set(node.pop("classes"))
        if "supplement" in classes and "main" not in classes:
            color = SYSTEM_STYLES[system]["supplement_color"]
        elif node["endpoint"]:
            color = SYSTEM_STYLES[system]["main_color"]
        else:
            color = SYSTEM_STYLES[system]["intermediate_color"]
        node["color"] = color
        node["radius"] = 1.05 if node["endpoint"] else 0.74
        nodes.append(node)
    nodes.sort(key=lambda item: int(item["canonical_resid"]))
    return routes, nodes, warnings


def pml_residue_label(node: dict[str, object]) -> str:
    return str(node["label"])


def html_doc(system: str, pdb_text: str, routes: list[dict[str, object]], nodes: list[dict[str, object]]) -> str:
    style = SYSTEM_STYLES[system]
    route_rows = []
    for route in routes:
        distance = route["mean_path_distance"]
        distance_text = "" if distance is None else f"; distance={distance:.3f}"
        support = f"{route['support_repeats']}/3 modal"
        route_rows.append(
            "<li>"
            f"<span class=\"swatch\" style=\"background:{html.escape(str(route['color']))}\"></span>"
            f"{html.escape(str(route['route_label']))}: "
            f"{html.escape(' -> '.join(route['path_labels']))} "
            f"({html.escape(str(route['pkg']))}; {support}{distance_text})"
            "</li>"
        )
    route_html = "\n".join(route_rows)
    return f"""<!doctype html>
<html>
<head>
  <meta charset="utf-8">
  <title>{html.escape(style["title"])}</title>
  <script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
  <style>
    body {{ margin: 0; font-family: Arial, sans-serif; color: #1e1e1e; }}
    #viewer {{ width: 100vw; height: 100vh; }}
    #legend {{
      position: absolute; left: 14px; top: 12px; z-index: 10;
      background: rgba(255,255,255,0.92); padding: 10px 12px;
      border: 1px solid #cfcfcf; border-radius: 4px; font-size: 13px;
      max-width: 560px; line-height: 1.35;
    }}
    #legend ul {{ margin: 7px 0 0 18px; padding: 0; }}
    #legend li {{ margin: 3px 0; }}
    .swatch {{ display:inline-block; width:12px; height:12px; margin-right:6px; vertical-align:middle; }}
    .note {{ margin-top: 7px; color: #444; font-size: 12px; }}
  </style>
</head>
<body>
  <div id="legend">
    <b>{html.escape(style["title"])}</b>
    <ul>
      {route_html}
    </ul>
    <div class="note">Routes are candidate communication paths from AlloViz-derived residue networks, not direct causal proof.</div>
  </div>
  <div id="viewer"></div>
  <script>
    const pdb = {js_string(pdb_text)};
    const routes = {json.dumps(routes)};
    const nodes = {json.dumps(nodes)};
    const viewer = $3Dmol.createViewer("viewer", {{backgroundColor: "white"}});
    viewer.addModel(pdb, "pdb");
    viewer.setStyle({{}}, {{cartoon: {{color: "lightgray", opacity: 0.68}}}});
    const coordByResid = {{}};
    for (const node of nodes) {{
      coordByResid[node.raw_resid] = node.coord;
      viewer.setStyle({{resi: node.raw_resid}}, {{
        sphere: {{radius: node.radius, color: node.color, opacity: 0.96}}
      }});
      viewer.addLabel(node.label, {{
        position: node.coord,
        fontSize: 13,
        fontColor: "black",
        backgroundColor: "white",
        backgroundOpacity: 0.68,
        borderThickness: 0.2,
        inFront: true
      }});
    }}
    for (const route of routes) {{
      for (const segment of route.segments) {{
        const start = coordByResid[segment.raw_a];
        const end = coordByResid[segment.raw_b];
        if (!start || !end) continue;
        viewer.addCylinder({{
          start: start,
          end: end,
          radius: segment.html_radius,
          color: route.color,
          opacity: 0.90,
          dashed: false
        }});
      }}
    }}
    viewer.zoomTo();
    viewer.rotate(20, {{x: 0, y: 1, z: 0}});
    viewer.render();
  </script>
</body>
</html>
"""


def pml_doc(system: str, routes: list[dict[str, object]], nodes: list[dict[str, object]]) -> str:
    style = SYSTEM_STYLES[system]
    lines = [
        f"# {style['title']}",
        "# Candidate communication routes inferred from AlloViz-derived networks.",
        "# Open in Windows PyMOL with File > Run Script.",
        "reinitialize",
        f'load "{SYSTEM_PDB_WINDOWS[system]}", receptor',
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
        f"set_color route_main, {style['pymol_main']}",
        f"set_color route_intermediate, {style['pymol_intermediate']}",
        f"set_color route_supplement, {style['pymol_supplement']}",
        "set_color route_label, [0.05, 0.05, 0.05]",
        "set label_color, route_label",
        "set label_size, 22",
        "set dash_gap, 0.00",
        "set dash_radius, 0.10",
        "set dash_round_ends, on",
        "# Line style controls:",
        "# - Each route segment below uses edge-score-derived dash_width/dash_radius.",
        "# - Larger dash_width/dash_radius means stronger mean edge score.",
        "# - set dash_gap, 0.00 gives a solid-looking line.",
        "",
        "# Route residue nodes.",
        "select active_route_nodes, none",
    ]
    for node in nodes:
        sel = safe_name(f"node_{node['label']}_{node['raw_resid']}")
        color = "route_supplement" if node["color"] == style["supplement_color"] else (
            "route_main" if node["endpoint"] else "route_intermediate"
        )
        scale = 0.72 if node["endpoint"] else 0.48
        lines.extend(
            [
                f"select {sel}, receptor and resi {node['raw_resid']}",
                f"show spheres, {sel} and name CA",
                f"set sphere_scale, {scale:.2f}, {sel} and name CA",
                f"color {color}, {sel} and name CA",
                f'label {sel} and name CA, "{pml_residue_label(node)}"',
                f"select active_route_nodes, active_route_nodes or {sel}",
            ]
        )
    lines.extend(["", "# Route segments as thick solid-looking CA-CA links."])
    for route_index, route in enumerate(routes, start=1):
        color = "route_supplement" if route["class_name"] == "supplement" else "route_main"
        label_path = list(route["path_labels"])
        lines.append(f"# {route['route_label']}: {' -> '.join(label_path)}")
        for segment in route["segments"]:
            segment_index = int(segment["index"])
            raw_a = int(segment["raw_a"])
            raw_b = int(segment["raw_b"])
            obj = safe_name(f"route_{route_index}_{segment_index}_{route['pair_label']}")
            width = float(segment["pymol_dash_width"])
            radius = float(segment["pymol_dash_radius"])
            gap = 0.00
            score = segment["strength_score"]
            score_text = "NA" if score is None else f"{float(score):.6g}"
            lines.extend(
                [
                    f"# segment {segment['label']}; method={segment['strength_pkg']}; mean_score={score_text}; support_repeats={segment['strength_support_repeats']}",
                    f"distance {obj}, receptor and resi {raw_a} and name CA, receptor and resi {raw_b} and name CA",
                    f"hide labels, {obj}",
                    f"color {color}, {obj}",
                    f"set dash_width, {width:.2f}, {obj}",
                    f"set dash_radius, {radius:.2f}, {obj}",
                    f"set dash_gap, {gap:.2f}, {obj}",
                ]
            )
    lines.extend(
        [
            "",
            "orient active_route_nodes",
            "zoom active_route_nodes, 12",
            "# Optional rendering commands:",
            "# ray 2400, 1800",
            f"# png {style['slug']}.png, dpi=300",
            "",
        ]
    )
    return "\n".join(lines)


def route_summary_rows(system: str, routes: list[dict[str, object]]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for route in routes:
        segments = list(route.get("segments", []))
        segment_strengths = []
        segment_widths = []
        for segment in segments:
            score = segment["strength_score"]
            score_text = "NA" if score is None else f"{float(score):.6g}"
            segment_strengths.append(f"{segment['label']}={score_text}")
            segment_widths.append(f"{segment['label']}={float(segment['pymol_dash_width']):.2f}")
        rows.append(
            {
                "system": system,
                "pkg": route["pkg"],
                "pair_label": route["pair_label"],
                "route_label": route["route_label"],
                "class_name": route["class_name"],
                "support_repeats": route["support_repeats"],
                "successful_repeats": route["successful_repeats"],
                "mean_path_distance": route["mean_path_distance"],
                "path_residues": " -> ".join(route["path_labels"]),
                "path_raw_resids": " -> ".join(str(item) for item in route["path_raw_resids"]),
                "segment_strengths": "; ".join(segment_strengths),
                "segment_dash_widths": "; ".join(segment_widths),
            }
        )
    return rows


def write_readme(manifest: dict[str, object]) -> None:
    lines = [
        "# C5aR1 Candidate Activation-Route Views",
        "",
        "This folder contains direct-view structure visualizations for candidate BM213 and C5apep communication routes.",
        "",
        "Use the HTML files for quick browser inspection. They embed the representative receptor PDB but load 3Dmol.js from `https://3Dmol.csb.pitt.edu`.",
        "",
        "Use the PML files in Windows PyMOL for manuscript-grade ray-traced images.",
        "",
        "| system | HTML | PyMOL |",
        "|---|---|---|",
    ]
    for item in manifest["systems"]:
        lines.append(f"| {item['system']} | `{Path(item['html']).name}` | `{Path(item['pml']).name}` |")
    lines.extend(
        [
            "",
            "Legend files: `activation_route_legend.png` and `activation_route_legend.svg`.",
            "",
            "Interpretation boundary: these are candidate communication routes inferred from AlloViz-derived residue networks and NetworkX shortest paths. They should be combined with microswitch/FES/structural evidence before being described as activation mechanisms.",
            "",
            "Route definitions are recorded in `activation_routes_summary.csv` and `activation_route_manifest.json`.",
        ]
    )
    (ACTIVE_DIR / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_legend_figure() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.linewidth": 0.8,
        }
    )
    fig, ax = plt.subplots(figsize=(5.8, 4.2), dpi=300)
    ax.set_xlim(0, 7)
    ax.set_ylim(0, 6.4)
    ax.axis("off")

    ax.text(
        0.3,
        5.95,
        "Activation-route colors",
        fontsize=15,
        fontweight="bold",
        ha="left",
        va="center",
    )
    ax.text(
        0.3,
        5.55,
        "Color encodes route class and residue role, not a continuous score.",
        fontsize=9.5,
        color="#444444",
        ha="left",
        va="center",
    )

    swatches = [
        ("BM213 main route", SYSTEM_STYLES["BM213"]["main_color"]),
        ("BM213 intermediate residue", SYSTEM_STYLES["BM213"]["intermediate_color"]),
        ("BM213 correlation supplement", SYSTEM_STYLES["BM213"]["supplement_color"]),
        ("C5apep main route", SYSTEM_STYLES["C5apep"]["main_color"]),
        ("C5apep intermediate residue", SYSTEM_STYLES["C5apep"]["intermediate_color"]),
        ("C5apep correlation supplement", SYSTEM_STYLES["C5apep"]["supplement_color"]),
    ]
    y = 4.85
    for label, color in swatches:
        ax.add_patch(Rectangle((0.38, y - 0.15), 0.40, 0.30, color=color, ec="#222222", lw=0.45))
        ax.text(0.98, y, label, ha="left", va="center", fontsize=10.5)
        y -= 0.62

    fig.tight_layout(pad=0.7)
    fig.savefig(ACTIVE_DIR / "activation_route_legend.png", dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(ACTIVE_DIR / "activation_route_legend.svg", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main() -> int:
    ACTIVE_DIR.mkdir(parents=True, exist_ok=True)
    consensus = pd.read_csv(TABLES / "networkx_consensus_routes.csv")
    mapping = pd.read_csv(TABLES / "residue_mapping_alloviz.csv")
    edge_scores = pd.read_csv(EDGE_SCORES)

    manifest: dict[str, object] = {
        "note": "Candidate activation-like route views from AlloViz-derived NetworkX consensus routes.",
        "inputs": {
            "routes": "tables/networkx_consensus_routes.csv",
            "mapping": "tables/residue_mapping_alloviz.csv",
            "pdbs": {key: str(path.relative_to(OUTROOT)) for key, path in SYSTEM_PDB.items()},
        },
        "systems": [],
        "warnings": [],
    }
    summary_rows: list[dict[str, object]] = []

    for system in ["BM213", "C5apep"]:
        pdb_path = SYSTEM_PDB[system]
        pdb_text = pdb_path.read_text(encoding="utf-8", errors="replace")
        coords = parse_ca_coords(pdb_path)
        routes, nodes, warnings = build_routes(system, consensus, mapping, coords)
        annotate_route_segments(routes, edge_scores)
        style = SYSTEM_STYLES[system]
        html_path = ACTIVE_DIR / f"{style['slug']}.html"
        pml_path = ACTIVE_DIR / f"{style['slug']}.pml"
        html_path.write_text(html_doc(system, pdb_text, routes, nodes), encoding="utf-8")
        pml_path.write_text(pml_doc(system, routes, nodes) + "\n", encoding="utf-8")
        summary_rows.extend(route_summary_rows(system, routes))
        manifest["systems"].append(
            {
                "system": system,
                "html": str(html_path.relative_to(OUTROOT)),
                "pml": str(pml_path.relative_to(OUTROOT)),
                "route_count": len(routes),
                "node_count": len(nodes),
            }
        )
        manifest["warnings"].extend(warnings)

    pd.DataFrame(summary_rows).to_csv(ACTIVE_DIR / "activation_routes_summary.csv", index=False)
    (ACTIVE_DIR / "activation_route_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    write_legend_figure()
    write_readme(manifest)
    print(f"Wrote activation route views under {ACTIVE_DIR}")
    for item in manifest["systems"]:
        print(f"{item['system']}: {item['html']} | {item['pml']}")
    if manifest["warnings"]:
        print("Warnings:")
        for warning in manifest["warnings"]:
            print(f"- {warning}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
