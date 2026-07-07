#!/usr/bin/env python3
from __future__ import annotations

import ast
import json
import math
import re
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import MDAnalysis as mda
import numpy as np
import pandas as pd
from MDAnalysis.lib.distances import capped_distance


ROOT = Path("/mnt/e/work/modeling/c5ar/md")
OUTROOT = ROOT / "article/analysis/alloviz"
RUNROOT = OUTROOT / "runs/formal"
TABLES = OUTROOT / "tables"
FIGURES = OUTROOT / "figures"
LOGS = OUTROOT / "logs"

THREE_TO_ONE = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "HSD": "H",
    "HSE": "H",
    "HSP": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}


@dataclass(frozen=True)
class SystemSpec:
    key: str
    pdb: str
    trajs: tuple[str, ...]
    receptor_sel: str
    ligand_sel: str | None
    role: str


SYSTEMS = [
    SystemSpec(
        key="apo",
        pdb="gpcr/apo/md-initial.pdb",
        trajs=(
            "gpcr/apo/2000ns.xtc",
            "gpcr/apo/second/2000ns.xtc",
            "gpcr/apo/third/2000ns.xtc",
        ),
        receptor_sel="protein",
        ligand_sel=None,
        role="basal reference",
    ),
    SystemSpec(
        key="BM213",
        pdb="gpcr/bm213/md-initial.pdb",
        trajs=(
            "gpcr/bm213/2wframe.xtc",
            "gpcr/bm213/second/2wframe.xtc",
            "gpcr/bm213/third/2wframe.xtc",
        ),
        receptor_sel="protein and not resid 282-289",
        ligand_sel="resid 282-289",
        role="G protein-biased / stronger agonist",
    ),
    SystemSpec(
        key="C5apep",
        pdb="gpcr/c5apep/md-initial.pdb",
        trajs=(
            "gpcr/c5apep/2wframe.xtc",
            "gpcr/c5apep/second/2wframe2.xtc",
            "gpcr/c5apep/third/2wframe3.xtc",
        ),
        receptor_sel="protein",
        ligand_sel="resid 284",
        role="partial efficacy",
    ),
]

KEY_RESIDUES = [
    "I91",
    "W102",
    "I116",
    "M120",
    "R134",
    "S171",
    "Y222",
    "W255",
    "F251",
    "Y258",
    "Y290",
    "N292",
    "N296",
    "Y300",
]

MOTIF_SHORT_BY_CANONICAL_RESID = {
    58: "IWI",
    69: "IWI",
    83: "IWI",
    87: "microswitch",
    101: "DRY",
    138: "pocket",
    189: "TM5-DRY",
    218: "PIF/CWxP",
    222: "CWxP",
    225: "CWxP",
    257: "TM7",
    259: "NPxxY",
    263: "NPxxY",
    267: "NPxxY",
}

def ensure_dirs() -> None:
    for path in (TABLES, FIGURES, LOGS):
        path.mkdir(parents=True, exist_ok=True)


def one_letter(resname: str) -> str:
    return THREE_TO_ONE.get(resname.upper(), "X")


def sequence_for(spec: SystemSpec) -> tuple[str, list[dict[str, object]]]:
    u = mda.Universe(str(ROOT / spec.pdb))
    residues = u.select_atoms(spec.receptor_sel).residues
    rows = []
    seq = []
    for idx, residue in enumerate(residues, start=1):
        aa = one_letter(residue.resname)
        seq.append(aa)
        rows.append(
            {
                "system": spec.key,
                "raw_resid": int(residue.resid),
                "resname": residue.resname,
                "aa": aa,
                "system_receptor_index": idx,
            }
        )
    return "".join(seq), rows


def nw_align(seq_a: str, seq_b: str) -> tuple[str, str]:
    match = 2
    mismatch = -1
    gap = -2
    n, m = len(seq_a), len(seq_b)
    score = np.zeros((n + 1, m + 1), dtype=float)
    trace = np.zeros((n + 1, m + 1), dtype=np.int8)
    for i in range(1, n + 1):
        score[i, 0] = i * gap
        trace[i, 0] = 1
    for j in range(1, m + 1):
        score[0, j] = j * gap
        trace[0, j] = 2
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diag = score[i - 1, j - 1] + (match if seq_a[i - 1] == seq_b[j - 1] else mismatch)
            up = score[i - 1, j] + gap
            left = score[i, j - 1] + gap
            best = max(diag, up, left)
            score[i, j] = best
            trace[i, j] = 0 if best == diag else (1 if best == up else 2)
    a_aln, b_aln = [], []
    i, j = n, m
    while i > 0 or j > 0:
        t = trace[i, j]
        if i > 0 and j > 0 and t == 0:
            a_aln.append(seq_a[i - 1])
            b_aln.append(seq_b[j - 1])
            i -= 1
            j -= 1
        elif i > 0 and (j == 0 or t == 1):
            a_aln.append(seq_a[i - 1])
            b_aln.append("-")
            i -= 1
        else:
            a_aln.append("-")
            b_aln.append(seq_b[j - 1])
            j -= 1
    return "".join(reversed(a_aln)), "".join(reversed(b_aln))


def load_key_mapping() -> pd.DataFrame:
    mapping_path = ROOT / "article/figure/tables/residue_mapping.csv"
    if not mapping_path.exists():
        return pd.DataFrame()
    df = pd.read_csv(mapping_path)
    df["trajectory_resi"] = df["trajectory_resi"].astype(str)
    return df


def build_residue_mapping() -> pd.DataFrame:
    seqs: dict[str, str] = {}
    residue_rows: dict[str, list[dict[str, object]]] = {}
    for spec in SYSTEMS:
        seq, rows = sequence_for(spec)
        seqs[spec.key] = seq
        residue_rows[spec.key] = rows

    ref_key = "BM213"
    ref_rows = residue_rows[ref_key]
    ref_idx_to_raw = {i + 1: int(row["raw_resid"]) for i, row in enumerate(ref_rows)}
    ref_raw_to_row = {int(row["raw_resid"]): row for row in ref_rows}

    raw_to_canonical: dict[tuple[str, int], int] = {}
    for row in ref_rows:
        raw_to_canonical[(ref_key, int(row["raw_resid"]))] = int(row["raw_resid"])

    for system in ("apo", "C5apep"):
        ref_aln, sys_aln = nw_align(seqs[ref_key], seqs[system])
        ref_pos = 0
        sys_pos = 0
        for ref_char, sys_char in zip(ref_aln, sys_aln):
            if ref_char != "-":
                ref_pos += 1
            if sys_char != "-":
                sys_pos += 1
            if ref_char == "-" or sys_char == "-":
                continue
            sys_raw = int(residue_rows[system][sys_pos - 1]["raw_resid"])
            raw_to_canonical[(system, sys_raw)] = ref_idx_to_raw[ref_pos]

    key_df = load_key_mapping()
    key_labels: dict[int, dict[str, str]] = {}
    if not key_df.empty:
        bm = key_df[(key_df["system"] == "BM213") & (key_df["panel_context"] != "ligand")]
        for _, row in bm.iterrows():
            traj_resi = str(row["trajectory_resi"])
            if not traj_resi.isdigit():
                continue
            canonical = int(traj_resi)
            key_labels[canonical] = {
                "biological_residue": str(row["biological_residue"]),
                "gpcrdb": str(row["bw_or_region"]),
                "key_label": str(row["label"]),
                "panel_context": str(row["panel_context"]),
            }
    manual_bm_key = {
        218: {"biological_residue": "F251", "gpcrdb": "6.44", "key_label": "F251", "panel_context": "CWxP/PIF"},
    }
    key_labels.update(manual_bm_key)

    rows: list[dict[str, object]] = []
    for system, sys_rows in residue_rows.items():
        for row in sys_rows:
            raw = int(row["raw_resid"])
            canonical = raw_to_canonical.get((system, raw))
            if canonical is None:
                continue
            ref_row = ref_raw_to_row.get(canonical, {})
            label_meta = key_labels.get(canonical, {})
            ref_resname = str(ref_row.get("resname", row["resname"]))
            default_label = f"{one_letter(ref_resname)}{canonical}"
            rows.append(
                {
                    "system": system,
                    "raw_resid": raw,
                    "resname": row["resname"],
                    "aa": row["aa"],
                    "system_receptor_index": row["system_receptor_index"],
                    "canonical_resid": canonical,
                    "canonical_resname": ref_resname,
                    "canonical_label": label_meta.get("key_label", default_label),
                    "biological_residue": label_meta.get("biological_residue", ""),
                    "gpcrdb": label_meta.get("gpcrdb", ""),
                    "panel_context": label_meta.get("panel_context", ""),
                    "is_key_residue": bool(label_meta),
                }
            )
    df = pd.DataFrame(rows)
    df.to_csv(TABLES / "residue_mapping_alloviz.csv", index=False)
    return df


def parse_numbers(value: object) -> list[int]:
    if pd.isna(value):
        return []
    text = str(value)
    return [int(x) for x in re.findall(r"(?<![A-Za-z])(\d+)(?![A-Za-z])", text)]


def parse_node_resid(value: object) -> int | None:
    nums = parse_numbers(value)
    return nums[-1] if nums else None


def parse_edge_resids(row: pd.Series) -> tuple[int | None, int | None]:
    unnamed = [col for col in row.index if str(col).startswith("Unnamed") or str(col).lower() in {"source", "target"}]
    if len(unnamed) >= 2:
        a = parse_node_resid(row[unnamed[0]])
        b = parse_node_resid(row[unnamed[1]])
        if a is not None and b is not None:
            return a, b

    first = row.iloc[0]
    if isinstance(first, str):
        try:
            parsed = ast.literal_eval(first)
            if isinstance(parsed, (tuple, list)) and len(parsed) == 2:
                return parse_node_resid(parsed[0]), parse_node_resid(parsed[1])
        except Exception:
            pass
    nums = parse_numbers(first)
    if len(nums) >= 2:
        return nums[-2], nums[-1]
    return None, None


def score_column(df: pd.DataFrame) -> str:
    if "btw" in df.columns:
        return "btw"
    numeric_cols = [col for col in df.columns if pd.api.types.is_numeric_dtype(df[col])]
    for col in numeric_cols:
        if not str(col).startswith("Unnamed"):
            return str(col)
    raise ValueError(f"No numeric score column in columns: {list(df.columns)}")


def load_run_exports(mapping: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
    map_lookup = {
        (row.system, int(row.raw_resid)): row
        for row in mapping.itertuples(index=False)
    }
    nodes: list[dict[str, object]] = []
    edges: list[dict[str, object]] = []
    warnings: list[str] = []

    for summary_path in sorted(RUNROOT.glob("*/*/*/run_summary.json")):
        system = summary_path.parents[2].name
        replicate = summary_path.parents[1].name
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        for export in summary.get("exports", []):
            pkg = export.get("pkg") or summary_path.parent.name
            element = export.get("element")
            path = Path(export["file"])
            if not path.exists():
                warnings.append(f"Missing export: {path}")
                continue
            df = pd.read_csv(path)
            try:
                score_col = score_column(df)
            except ValueError as exc:
                warnings.append(str(exc))
                continue
            if element == "nodes":
                id_col = df.columns[0]
                for _, row in df.iterrows():
                    raw = parse_node_resid(row[id_col])
                    if raw is None:
                        continue
                    mapped = map_lookup.get((system, raw))
                    if mapped is None:
                        continue
                    nodes.append(
                        {
                            "system": system,
                            "replicate": replicate,
                            "pkg": pkg,
                            "filtering": export.get("filtering", "Spatially_distant"),
                            "metric": score_col,
                            "raw_resid": raw,
                            "canonical_resid": int(mapped.canonical_resid),
                            "canonical_label": mapped.canonical_label,
                            "biological_residue": mapped.biological_residue,
                            "gpcrdb": mapped.gpcrdb,
                            "panel_context": mapped.panel_context,
                            "score": float(row[score_col]),
                        }
                    )
            elif element == "edges":
                for _, row in df.iterrows():
                    raw_a, raw_b = parse_edge_resids(row)
                    if raw_a is None or raw_b is None:
                        continue
                    map_a = map_lookup.get((system, raw_a))
                    map_b = map_lookup.get((system, raw_b))
                    if map_a is None or map_b is None:
                        continue
                    ca, cb = int(map_a.canonical_resid), int(map_b.canonical_resid)
                    if ca == cb:
                        continue
                    if ca > cb:
                        ca, cb = cb, ca
                        map_a, map_b = map_b, map_a
                    edges.append(
                        {
                            "system": system,
                            "replicate": replicate,
                            "pkg": pkg,
                            "filtering": export.get("filtering", "Spatially_distant"),
                            "metric": score_col,
                            "raw_resid_a": raw_a,
                            "raw_resid_b": raw_b,
                            "canonical_resid_a": ca,
                            "canonical_resid_b": cb,
                            "canonical_label_a": map_a.canonical_label,
                            "canonical_label_b": map_b.canonical_label,
                            "edge_label": f"{map_a.canonical_label}-{map_b.canonical_label}",
                            "score": float(row[score_col]),
                        }
                    )

    node_df = pd.DataFrame(nodes)
    edge_df = pd.DataFrame(edges)
    if not node_df.empty:
        node_df.to_csv(TABLES / "node_scores_long.csv", index=False)
    if not edge_df.empty:
        edge_df.to_csv(TABLES / "edge_scores_long.csv", index=False)
    return node_df, edge_df, warnings


def write_delta_tables(node_df: pd.DataFrame, edge_df: pd.DataFrame) -> dict[str, pd.DataFrame]:
    outputs: dict[str, pd.DataFrame] = {}
    comparisons = [("BM213", "apo"), ("C5apep", "apo"), ("BM213", "C5apep")]

    for a, b in comparisons:
        if not node_df.empty:
            pivot = node_df.pivot_table(
                index=["pkg", "canonical_resid", "canonical_label", "biological_residue", "gpcrdb", "panel_context"],
                columns="system",
                values="score",
                aggfunc="mean",
            ).reset_index()
            if a in pivot.columns and b in pivot.columns:
                out = pivot.dropna(subset=[a, b]).copy()
                out[f"delta_{a}_minus_{b}"] = out[a] - out[b]
                out = out.sort_values(f"delta_{a}_minus_{b}", key=lambda s: s.abs(), ascending=False)
                filename = f"delta_nodes_{a}_vs_{b}.csv"
                out.to_csv(TABLES / filename, index=False)
                outputs[filename] = out

        if not edge_df.empty:
            pivot = edge_df.pivot_table(
                index=["pkg", "canonical_resid_a", "canonical_resid_b", "edge_label"],
                columns="system",
                values="score",
                aggfunc="mean",
            ).reset_index()
            if a in pivot.columns and b in pivot.columns:
                out = pivot.dropna(subset=[a, b]).copy()
                out[f"delta_{a}_minus_{b}"] = out[a] - out[b]
                out = out.sort_values(f"delta_{a}_minus_{b}", key=lambda s: s.abs(), ascending=False)
                filename = f"delta_edges_{a}_vs_{b}.csv"
                out.to_csv(TABLES / filename, index=False)
                outputs[filename] = out
    return outputs


def ligand_contacts(stride: int = 10, cutoff: float = 4.5) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    mapping = pd.read_csv(TABLES / "residue_mapping_alloviz.csv")
    map_lookup = {(row.system, int(row.raw_resid)): row for row in mapping.itertuples(index=False)}

    for spec in SYSTEMS:
        if not spec.ligand_sel:
            continue
        for rep_idx, traj in enumerate(spec.trajs, start=1):
            u = mda.Universe(str(ROOT / spec.pdb), str(ROOT / traj))
            receptor = u.select_atoms(spec.receptor_sel)
            ligand = u.select_atoms(spec.ligand_sel)
            if ligand.n_atoms == 0:
                continue
            counts: defaultdict[int, int] = defaultdict(int)
            total = 0
            for ts in u.trajectory[::stride]:
                total += 1
                pairs = capped_distance(
                    ligand.positions,
                    receptor.positions,
                    max_cutoff=cutoff,
                    box=ts.dimensions,
                    return_distances=False,
                )
                if len(pairs) == 0:
                    continue
                receptor_atom_ids = np.unique(pairs[:, 1])
                resids = set(int(atom.residue.resid) for atom in receptor[receptor_atom_ids])
                for resid in resids:
                    counts[resid] += 1
            for resid, count in counts.items():
                mapped = map_lookup.get((spec.key, resid))
                if mapped is None:
                    continue
                rows.append(
                    {
                        "system": spec.key,
                        "replicate": f"rep{rep_idx}",
                        "trajectory": traj,
                        "ligand_selection": spec.ligand_sel,
                        "cutoff_A": cutoff,
                        "stride": stride,
                        "sampled_frames": total,
                        "raw_resid": resid,
                        "canonical_resid": int(mapped.canonical_resid),
                        "canonical_label": mapped.canonical_label,
                        "biological_residue": mapped.biological_residue,
                        "gpcrdb": mapped.gpcrdb,
                        "contact_frames": count,
                        "contact_frequency": count / total if total else math.nan,
                    }
                )

    df = pd.DataFrame(rows)
    if not df.empty:
        df.to_csv(TABLES / "ligand_receptor_contact_summary.csv", index=False)
        agg = (
            df.groupby(["system", "canonical_resid", "canonical_label", "biological_residue", "gpcrdb"], dropna=False)
            .agg(mean_contact_frequency=("contact_frequency", "mean"), max_contact_frequency=("contact_frequency", "max"))
            .reset_index()
            .sort_values(["system", "mean_contact_frequency"], ascending=[True, False])
        )
        agg.to_csv(TABLES / "ligand_receptor_contact_summary_by_residue.csv", index=False)
    return df


def key_residue_summary(node_df: pd.DataFrame, ligand_df: pd.DataFrame) -> pd.DataFrame:
    if node_df.empty:
        return pd.DataFrame()
    key = node_df[node_df["canonical_label"].isin(KEY_RESIDUES) | node_df["biological_residue"].isin(KEY_RESIDUES)].copy()
    if key.empty:
        key = node_df[node_df["panel_context"].astype(str) != ""].copy()
    pivot = key.pivot_table(
        index=["pkg", "canonical_resid", "canonical_label", "biological_residue", "gpcrdb", "panel_context"],
        columns="system",
        values="score",
        aggfunc="mean",
    ).reset_index()
    for a, b in [("BM213", "apo"), ("C5apep", "apo"), ("BM213", "C5apep")]:
        if a in pivot.columns and b in pivot.columns:
            pivot[f"delta_{a}_minus_{b}"] = pivot[a] - pivot[b]
    if not ligand_df.empty:
        lig = (
            ligand_df.groupby(["system", "canonical_resid"])["contact_frequency"]
            .mean()
            .reset_index()
            .pivot(index="canonical_resid", columns="system", values="contact_frequency")
            .add_prefix("ligand_contact_")
            .reset_index()
        )
        pivot = pivot.merge(lig, on="canonical_resid", how="left")
    pivot.to_csv(TABLES / "key_residue_network_summary.csv", index=False)
    return pivot


def clean_label(value: object) -> str:
    if pd.isna(value):
        return ""
    text = str(value).strip()
    if text.lower() == "nan":
        return ""
    return text


def residue_display_label(row: pd.Series) -> str:
    residue = clean_label(row.get("biological_residue")) or clean_label(row.get("canonical_label"))
    if not residue:
        residue = str(int(row["canonical_resid"]))
    gpcrdb = clean_label(row.get("gpcrdb"))
    if gpcrdb:
        return f"{residue} ({gpcrdb})"
    return residue


def motif_short_label(row: pd.Series) -> str:
    canonical_resid = int(row["canonical_resid"])
    motif = MOTIF_SHORT_BY_CANONICAL_RESID.get(canonical_resid, "")
    if motif:
        return motif
    context = clean_label(row.get("panel_context"))
    if not context:
        return ""
    return context.split("/")[0].strip()


def heatmap_axis_label(row: pd.Series) -> str:
    residue = residue_display_label(row)
    motif = motif_short_label(row)
    if not motif:
        return ""
    return f"{residue}-{motif}"


def prepare_motif_heatmap_rows(key_summary: pd.DataFrame, preferred_pkg: str) -> pd.DataFrame:
    df = key_summary[key_summary["pkg"] == preferred_pkg].copy()
    df["motif_axis_label"] = df.apply(heatmap_axis_label, axis=1)
    return df[df["motif_axis_label"] != ""].copy()


def plot_annotated_heatmap(
    df: pd.DataFrame,
    value_cols: list[str],
    xlabels: list[str],
    title: str,
    colorbar_label: str,
    filename: str,
    cmap_name: str,
    center_zero: bool,
) -> None:
    if df.empty or not value_cols:
        return
    labels = df["motif_axis_label"].tolist()
    matrix = df[value_cols].to_numpy(dtype=float)
    masked = np.ma.masked_invalid(matrix)
    if center_zero:
        vmax = np.nanmax(np.abs(matrix)) if matrix.size else 1.0
        vmax = max(vmax, 1e-6)
        vmin = -vmax
    else:
        vmin = 0.0
        vmax = np.nanmax(matrix) if matrix.size else 1.0
        vmax = max(vmax, 1e-6)

    cmap = plt.get_cmap(cmap_name).copy()
    cmap.set_bad("#E6E6E6")

    height = max(6.4, 0.62 * len(labels) + 1.3)
    fig = plt.figure(figsize=(15.4, height))
    gs = fig.add_gridspec(1, 2, width_ratios=[4.2, 0.22], wspace=0.08)
    ax = fig.add_subplot(gs[0, 0])
    cax = fig.add_subplot(gs[0, 1])

    im = ax.imshow(masked, cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto")
    ax.set_yticks(np.arange(len(labels)), labels=labels)
    ax.set_xticks(np.arange(len(xlabels)), labels=xlabels, rotation=0, ha="center")
    ax.tick_params(axis="y", labelsize=20)
    ax.tick_params(axis="x", labelsize=26, pad=10)
    ax.set_title(title, fontsize=24, pad=18)
    ax.set_xticks(np.arange(-0.5, len(xlabels), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(labels), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.7)
    ax.tick_params(which="minor", bottom=False, left=False)

    cbar = fig.colorbar(im, cax=cax)
    cbar.set_label(colorbar_label, fontsize=22)
    cbar.ax.tick_params(labelsize=18)
    fig.savefig(FIGURES / filename, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_node_delta_heatmap(key_summary: pd.DataFrame) -> None:
    if key_summary.empty:
        return
    pkg_order = list(dict.fromkeys(key_summary["pkg"].tolist()))
    preferred_pkg = "GetContacts" if "GetContacts" in pkg_order else pkg_order[0]
    df = prepare_motif_heatmap_rows(key_summary, preferred_pkg)
    delta_cols = [col for col in df.columns if col.startswith("delta_")]
    if not delta_cols:
        return
    plot_annotated_heatmap(
        df=df,
        value_cols=delta_cols,
        xlabels=[c.replace("delta_", "").replace("_minus_", "-") for c in delta_cols],
        title=f"AlloViz node betweenness deltas ({preferred_pkg})",
        colorbar_label="delta btw",
        filename="alloviz_node_delta_heatmap.png",
        cmap_name="coolwarm",
        center_zero=True,
    )


def plot_node_absolute_heatmap(key_summary: pd.DataFrame) -> None:
    if key_summary.empty:
        return
    pkg_order = list(dict.fromkeys(key_summary["pkg"].tolist()))
    preferred_pkg = "GetContacts" if "GetContacts" in pkg_order else pkg_order[0]
    df = prepare_motif_heatmap_rows(key_summary, preferred_pkg)
    state_cols = [col for col in ("apo", "BM213", "C5apep") if col in df.columns]
    if not state_cols:
        return
    plot_annotated_heatmap(
        df=df,
        value_cols=state_cols,
        xlabels=state_cols,
        title=f"AlloViz node betweenness scores ({preferred_pkg})",
        colorbar_label="btw score",
        filename="alloviz_node_absolute_heatmap.png",
        cmap_name="viridis",
        center_zero=False,
    )


def plot_top_edges(edge_df: pd.DataFrame) -> None:
    if edge_df.empty:
        return
    df = edge_df[edge_df["pkg"] == "GetContacts"].copy()
    if df.empty:
        df = edge_df.copy()
    rows = []
    for system in ("apo", "BM213", "C5apep"):
        sub = df[df["system"] == system].nlargest(8, "score")
        for _, row in sub.iterrows():
            rows.append({"system": system, "edge": row["edge_label"], "score": row["score"]})
    plot_df = pd.DataFrame(rows)
    if plot_df.empty:
        return
    fig, axes = plt.subplots(1, 3, figsize=(14, 5), sharex=False)
    for ax, system in zip(axes, ("apo", "BM213", "C5apep")):
        sub = plot_df[plot_df["system"] == system].sort_values("score")
        ax.barh(sub["edge"], sub["score"], color="#4C78A8")
        ax.set_title(system)
        ax.set_xlabel("btw")
    fig.suptitle("Top AlloViz edge scores by state")
    fig.tight_layout()
    fig.savefig(FIGURES / "alloviz_top_edges_by_state.png", dpi=300)
    plt.close(fig)


def plot_key_residue_network(edge_df: pd.DataFrame, mapping: pd.DataFrame) -> None:
    if edge_df.empty:
        return
    key_map = mapping[mapping["is_key_residue"]].copy()
    key_resids = set(key_map["canonical_resid"].astype(int))
    df = edge_df[(edge_df["pkg"] == "GetContacts") & (edge_df["canonical_resid_a"].isin(key_resids) | edge_df["canonical_resid_b"].isin(key_resids))].copy()
    if df.empty:
        df = edge_df[edge_df["pkg"] == "GetContacts"].nlargest(30, "score").copy()
    df = df.nlargest(40, "score")
    nodes = sorted(set(df["canonical_label_a"]).union(set(df["canonical_label_b"])))
    if not nodes:
        return
    angles = np.linspace(0, 2 * np.pi, len(nodes), endpoint=False)
    pos = {node: np.array([np.cos(a), np.sin(a)]) for node, a in zip(nodes, angles)}
    color_map = {"apo": "#888888", "BM213": "#D55E00", "C5apep": "#0072B2"}
    key_labels = set(key_map["canonical_label"])
    fig, ax = plt.subplots(figsize=(21.0, 17.0))
    for _, row in df.iterrows():
        a, b = row["canonical_label_a"], row["canonical_label_b"]
        pa, pb = pos[a], pos[b]
        ax.plot(
            [pa[0], pb[0]],
            [pa[1], pb[1]],
            color=color_map.get(row["system"], "#444444"),
            alpha=0.18 + 0.52 * min(float(row["score"]) / max(df["score"].max(), 1e-9), 1),
            linewidth=0.8 + 3.0 * min(float(row["score"]) / max(df["score"].max(), 1e-9), 1),
            zorder=1,
        )

    for node, p in pos.items():
        is_key = node in key_labels
        ax.scatter(
            p[0],
            p[1],
            s=390 if is_key else 165,
            color="#FFFFFF" if is_key else "#F3F3F3",
            edgecolor="#222222" if is_key else "#555555",
            linewidth=1.6 if is_key else 1.1,
            zorder=3,
        )

    for node, p in pos.items():
        label_pos = p * 1.28
        x, y = float(label_pos[0]), float(label_pos[1])
        if abs(float(p[0])) < 0.2:
            ha = "center"
        else:
            ha = "left" if p[0] > 0 else "right"
        if abs(float(p[1])) < 0.2:
            va = "center"
        else:
            va = "bottom" if p[1] > 0 else "top"
        ax.text(
            x,
            y,
            node,
            ha=ha,
            va=va,
            fontsize=34,
            fontweight="normal",
            bbox={"boxstyle": "round,pad=0.12", "facecolor": "white", "edgecolor": "none", "alpha": 0.82},
            zorder=4,
        )

    for system, color in color_map.items():
        ax.plot([], [], color=color, label=system, linewidth=4)
    ax.set_title("AlloViz key-residue network", fontsize=40, pad=22)
    ax.legend(frameon=False, loc="lower center", bbox_to_anchor=(0.5, -0.045), ncol=3, fontsize=36, handlelength=1.7, columnspacing=1.6)
    ax.set_xlim(-1.9, 1.9)
    ax.set_ylim(-1.95, 1.95)
    ax.set_axis_off()
    ax.set_facecolor("white")
    fig.tight_layout()
    fig.savefig(FIGURES / "alloviz_key_residue_network.png", dpi=600, bbox_inches="tight")
    fig.savefig(FIGURES / "alloviz_key_residue_network.pdf", bbox_inches="tight")
    plt.close(fig)


def write_report(node_df: pd.DataFrame, edge_df: pd.DataFrame, ligand_df: pd.DataFrame, warnings: list[str]) -> None:
    status_path = LOGS / "run_status_summary.json"
    statuses = json.loads(status_path.read_text(encoding="utf-8")) if status_path.exists() else []
    successes = [s for s in statuses if s.get("status") in {"success", "skipped_existing"} and s.get("phase") == "formal"]
    failures = [s for s in statuses if s.get("status") == "failed"]

    lines = [
        "# AlloViz C5aR1 Network Analysis",
        "",
        f"- Generated: `{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}`",
        "- Scope: apo, BM213, and C5apep receptor-only C5aR1 networks plus separate ligand contact layer.",
        "- Interpretation: dynamic network evidence and candidate allosteric communication routes, not direct causal proof.",
        "",
        "## Inputs",
        "",
        "- Trajectory source: `gpcr/dry_independent_replicates.md`, with C5apep rep3 overridden to `gpcr/c5apep/third/2wframe3.xtc`.",
        "- Receptor selections: apo `protein`; BM213 `protein and not resid 282-289`; C5apep `protein`.",
        "- Ligand layer: BM213 `resid 282-289`; C5apep `resid 284`.",
        "",
        "## AlloViz Runs",
        "",
        f"- Formal successful/skipped-existing runs: `{len(successes)}`",
        f"- Failed runs recorded: `{len(failures)}`",
    ]
    if successes:
        lines.extend(["", "| system | replicate | package | status |", "|---|---|---|---|"])
        for item in successes:
            lines.append(
                f"| {item.get('system')} | {item.get('replicate', 'NA')} | "
                f"{item.get('pkg')} | {item.get('status')} |"
            )
    if failures:
        lines.extend(["", "See `logs/method_failures.md` for failures."])

    lines.extend(
        [
            "",
            "## Output Tables",
            "",
            "- `tables/residue_mapping_alloviz.csv`: raw residue to canonical BM213-aligned receptor mapping.",
            "- `tables/node_scores_long.csv` and `tables/edge_scores_long.csv`: long-format AlloViz node/edge betweenness scores.",
            "- `tables/delta_nodes_*.csv` and `tables/delta_edges_*.csv`: common-residue/network deltas for BM213 vs apo, C5apep vs apo, and BM213 vs C5apep.",
            "- `tables/key_residue_network_summary.csv`: pathway-focused residues including Y258/S171/M120/W255/F251/R134/Y222/NPxxY/TM7 markers.",
            "- `tables/ligand_receptor_contact_summary*.csv`: ligand input-layer contact frequencies.",
            "",
            "## Figures",
            "",
            "- `figures/alloviz_node_delta_heatmap.png`",
            "- `figures/alloviz_node_absolute_heatmap.png`",
            "- `figures/alloviz_key_residue_network.png`",
            "- `figures/alloviz_top_edges_by_state.png`",
            "",
            "## Official AlloViz View Outputs",
            "",
            "- Use `scripts/generate_official_views.py` for tutorial-style `system.view(...)` and `AlloViz.Delta(...).view(...)` outputs.",
            "- Outputs are written under `official_views/` and include status JSON files plus `delta_alignment.txt` when Delta alignment succeeds.",
            "- These structure views are intended for visual panels; manuscript statistics should still come from the repeat-aware CSV tables.",
            "",
            "## Notes For Manuscript Use",
            "",
            "- Use wording such as `ligand-dependent network reorganization` or `candidate communication routes`.",
            "- Do not claim AlloViz proves a causal activation path by itself.",
            "- Combine these tables with existing microswitch/FES evidence before making claims about partial versus stronger agonism.",
        ]
    )
    if warnings:
        lines.extend(["", "## Parser Warnings", ""])
        for warning in warnings:
            lines.append(f"- {warning}")

    if not node_df.empty:
        top_nodes = node_df[node_df["pkg"] == "GetContacts"].sort_values("score", ascending=False).head(10)
        if not top_nodes.empty:
            lines.extend(["", "## Top GetContacts Nodes", "", "| system | residue | score |", "|---|---|---:|"])
            for _, row in top_nodes.iterrows():
                lines.append(f"| {row['system']} | {row['canonical_label']} | {row['score']:.4g} |")
    if not ligand_df.empty:
        lig = (
            ligand_df.groupby(["system", "canonical_label"], dropna=False)["contact_frequency"]
            .mean()
            .reset_index()
            .sort_values(["system", "contact_frequency"], ascending=[True, False])
            .groupby("system")
            .head(8)
        )
        lines.extend(["", "## Top Ligand Contacts", "", "| system | residue | mean frequency |", "|---|---|---:|"])
        for _, row in lig.iterrows():
            lines.append(f"| {row['system']} | {row['canonical_label']} | {row['contact_frequency']:.3f} |")

    (OUTROOT / "report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    ensure_dirs()
    mapping = build_residue_mapping()
    node_df, edge_df, warnings = load_run_exports(mapping)
    if node_df.empty and edge_df.empty:
        raise SystemExit("No formal AlloViz exports found. Run run_alloviz_methods.py --phase formal first.")
    write_delta_tables(node_df, edge_df)
    ligand_df = ligand_contacts()
    key_summary = key_residue_summary(node_df, ligand_df)
    plot_node_delta_heatmap(key_summary)
    plot_node_absolute_heatmap(key_summary)
    plot_top_edges(edge_df)
    plot_key_residue_network(edge_df, mapping)
    if warnings:
        (LOGS / "postprocess_warnings.txt").write_text("\n".join(warnings) + "\n", encoding="utf-8")
    write_report(node_df, edge_df, ligand_df, warnings)
    print(f"Wrote outputs under {OUTROOT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
