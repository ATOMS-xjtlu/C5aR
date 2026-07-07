#!/usr/bin/env python3
"""Summarize activation-core residue changes along BM213/C5apep SOM transitions."""

from __future__ import annotations

import math
import os
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-c5ar-sommd")

SYSTEMS = ("BM213", "C5apep")
EDGE_THRESHOLD = 5
CORE_RESIDUES = [
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

PALETTE = {
    "blue_main": "#0F4D92",
    "blue_secondary": "#3775BA",
    "neutral_light": "#CFCECE",
    "neutral_mid": "#767676",
    "neutral_dark": "#4D4D4D",
    "neutral_black": "#272727",
    "cyan": "#22D7E6",
    "warm_orange": "#E28E2C",
    "pastel_soft": "#B4C0E4",
    "pastel_bg": "#E0E0F0",
}


def label_to_nominal_number(label: str) -> int | None:
    match = re.search(r"(\d+)$", label)
    if not match:
        return None
    return int(match.group(1))


def label_to_canonical_resid(label: str) -> int | None:
    """Map display labels such as W255 to the unified C5aR residue id.

    The formal PyMOL structures were sequence-unified to common receptor residue
    ids. For the activation labels used here, common id = display number - 33.
    The core labels are also checked against activation_residue_mapping.csv.
    """

    number = label_to_nominal_number(label)
    if number is None:
        return None
    return number - 33


def parse_route_residues(route_summary: pd.DataFrame, system: str) -> set[str]:
    labels: set[str] = set()
    if route_summary.empty:
        return labels
    subset = route_summary[route_summary["system"] == system]
    for path_text in subset["path_residues"].dropna():
        for token in str(path_text).split("->"):
            token = token.strip()
            if token:
                labels.add(token)
    return labels


def find_representative_pdb(formal_dir: Path, system: str, cluster: int) -> tuple[Path | None, str]:
    archive_dir = formal_dir / "archive_non_dominant_cluster_pdbs"
    neuron_dir = formal_dir / "neuron_representatives"

    search_sets = [
        ("cluster_level_system_specific", sorted(formal_dir.glob(f"cluster_{cluster:02d}_{system}_*.pdb"))),
        ("archive_non_dominant_cluster", sorted(archive_dir.glob(f"cluster_{cluster:02d}_{system}_*.pdb"))),
        ("neuron_system_specific", sorted(neuron_dir.glob(f"neuron_*_C{cluster}_{system}_*.pdb"))),
        ("fallback_cluster_other_system", sorted(formal_dir.glob(f"cluster_{cluster:02d}_*.pdb"))),
    ]
    for source, paths in search_sets:
        if paths:
            return paths[0], source
    return None, "missing"


def read_ca_coords(pdb_path: Path) -> dict[int, tuple[str, np.ndarray]]:
    coords: dict[int, tuple[str, np.ndarray]] = {}
    with pdb_path.open() as handle:
        for line in handle:
            if not line.startswith("ATOM"):
                continue
            if line[12:16].strip() != "CA":
                continue
            try:
                resid = int(line[22:26])
                xyz = np.array(
                    [float(line[30:38]), float(line[38:46]), float(line[46:54])],
                    dtype=float,
                )
            except ValueError:
                continue
            resn = line[17:20].strip()
            coords[resid] = (resn, xyz)
    return coords


def kabsch_align(mobile: np.ndarray, target: np.ndarray) -> np.ndarray:
    mob_mean = mobile.mean(axis=0)
    tar_mean = target.mean(axis=0)
    mob_centered = mobile - mob_mean
    tar_centered = target - tar_mean
    covariance = mob_centered.T @ tar_centered
    v, _, wt = np.linalg.svd(covariance)
    correction = np.eye(3)
    correction[2, 2] = np.sign(np.linalg.det(v @ wt))
    rotation = v @ correction @ wt
    return mob_centered @ rotation + tar_mean


def safe_float(value: float) -> str:
    if pd.isna(value):
        return "NA"
    return f"{value:.2f}"


def apply_publication_style(font_size: float = 7.0) -> None:
    import matplotlib.pyplot as plt

    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
    plt.rcParams["svg.fonttype"] = "none"
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["font.size"] = font_size
    plt.rcParams["axes.linewidth"] = 0.7
    plt.rcParams["axes.spines.right"] = False
    plt.rcParams["axes.spines.top"] = False
    plt.rcParams["xtick.major.width"] = 0.7
    plt.rcParams["ytick.major.width"] = 0.7
    plt.rcParams["xtick.major.size"] = 2.5
    plt.rcParams["ytick.major.size"] = 2.5
    plt.rcParams["legend.frameon"] = False


def make_heatmap(residue_table: pd.DataFrame, out_png: Path, show_panel_labels: bool = True) -> None:
    import matplotlib as mpl
    import matplotlib.pyplot as plt

    apply_publication_style(font_size=8.2)
    cmap = mpl.colors.LinearSegmentedColormap.from_list(
        "nature_sequential_blue",
        ["#F8F8FA", PALETTE["pastel_bg"], PALETTE["pastel_soft"], PALETTE["blue_secondary"], PALETTE["blue_main"]],
    )
    cmap.set_bad("#FFFFFF")

    matrices: dict[str, pd.DataFrame] = {}
    column_labels: dict[str, list[str]] = {}
    for system in SYSTEMS:
        subset = residue_table[residue_table["system"] == system].copy()
        transitions = (
            subset[["from_cluster", "to_cluster", "transitions"]]
            .drop_duplicates()
            .sort_values(["transitions", "from_cluster", "to_cluster"], ascending=[False, True, True])
        )
        transition_keys = [(int(row.from_cluster), int(row.to_cluster)) for row in transitions.itertuples(index=False)]
        if system == "BM213" and (5, 1) in transition_keys and (1, 5) in transition_keys:
            transition_keys.remove((5, 1))
            transition_keys.insert(transition_keys.index((1, 5)), (5, 1))
        label_map = {
            (int(row.from_cluster), int(row.to_cluster)): f"C{int(row.from_cluster)}->C{int(row.to_cluster)}"
            for row in transitions.itertuples(index=False)
        }
        ordered_labels = [label_map[key] for key in transition_keys]
        subset["transition_label"] = [
            label_map[(int(row.from_cluster), int(row.to_cluster))] for row in subset.itertuples(index=False)
        ]
        matrix = subset.pivot_table(
            index="residue",
            columns="transition_label",
            values="mean_abs_delta_angstrom",
            aggfunc="first",
        )
        matrices[system] = matrix.reindex(CORE_RESIDUES)[ordered_labels]
        column_labels[system] = ordered_labels

    vmax = max(float(np.nanmax(matrix.values)) for matrix in matrices.values() if matrix.size)
    vmax = math.ceil(vmax * 10) / 10
    norm = mpl.colors.Normalize(vmin=0, vmax=vmax)

    fig = plt.figure(figsize=(7.75, 3.85), constrained_layout=True)
    grid = fig.add_gridspec(1, 3, width_ratios=[8, 5, 0.28], wspace=0.18)
    axes = [fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1])]
    cax = fig.add_subplot(grid[0, 2])

    panel_colors = {"BM213": PALETTE["neutral_black"], "C5apep": PALETTE["neutral_black"]}
    last_image = None
    for panel_label, ax, system in zip(("a", "b"), axes, SYSTEMS):
        matrix = matrices[system]
        if matrix.shape[1] == 0:
            ax.axis("off")
            ax.set_title(system)
            continue
        last_image = ax.imshow(matrix.values, aspect="auto", cmap=cmap, norm=norm, interpolation="nearest")
        ax.set_title(system, color=panel_colors[system], fontsize=10.2, fontweight="bold", pad=5)
        ax.set_xticks(np.arange(matrix.shape[1]))
        ax.set_xticklabels(matrix.columns, rotation=0, ha="center", fontsize=7.0)
        ax.set_yticks(np.arange(matrix.shape[0]))
        ax.set_yticklabels(matrix.index, fontsize=8.1)
        ax.set_xlabel("Dominant cluster transition", fontsize=8.4, labelpad=4)
        if system == "BM213":
            ax.set_ylabel("Activation-core residue", fontsize=8.4, labelpad=4)
        else:
            ax.set_ylabel("")
        ax.set_xticks(np.arange(-0.5, matrix.shape[1], 1), minor=True)
        ax.set_yticks(np.arange(-0.5, matrix.shape[0], 1), minor=True)
        ax.grid(which="minor", color="#FFFFFF", linewidth=0.6)
        ax.tick_params(which="minor", bottom=False, left=False)
        ax.tick_params(axis="both", colors=PALETTE["neutral_black"], pad=2)
        for spine in ax.spines.values():
            spine.set_color(PALETTE["neutral_dark"])
            spine.set_linewidth(0.7)
        if show_panel_labels:
            ax.text(
                -0.12,
                1.03,
                panel_label,
                transform=ax.transAxes,
                ha="left",
                va="bottom",
                fontsize=9.8,
                fontweight="bold",
                color=PALETTE["neutral_black"],
            )
    if last_image is not None:
        cbar = fig.colorbar(last_image, cax=cax)
        cbar.set_label(r"Mean $|\Delta d|$ per pair ($\mathrm{\AA}$)", fontsize=8.3, labelpad=5)
        cbar.ax.tick_params(labelsize=7.8, width=0.7, length=2.5)
        cbar.outline.set_linewidth(0.7)
        cbar.outline.set_edgecolor(PALETTE["neutral_dark"])
    fig.savefig(out_png, dpi=600, bbox_inches="tight", facecolor="white")
    fig.savefig(out_png.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    fig.savefig(out_png.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


def make_summary_barplot(system_residue_summary: pd.DataFrame, out_png: Path) -> None:
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4), constrained_layout=True)
    colors = {"BM213": "#1f77b4", "C5apep": "#ff8c00"}
    for ax, system in zip(axes, SYSTEMS):
        subset = (
            system_residue_summary[system_residue_summary["system"] == system]
            .sort_values("weighted_mean_abs_delta_angstrom", ascending=False)
            .head(10)
        )
        ax.barh(subset["residue"], subset["weighted_mean_abs_delta_angstrom"], color=colors[system])
        ax.invert_yaxis()
        ax.set_title(system)
        ax.set_xlabel("transition-weighted mean |distance change| (A)")
        ax.set_ylabel("")
    fig.savefig(out_png, dpi=300)
    plt.close(fig)


def write_pml(
    system: str,
    edges: pd.DataFrame,
    inventory: pd.DataFrame,
    system_residue_summary: pd.DataFrame,
    residue_to_canonical: dict[str, int],
    out_pml: Path,
) -> None:
    top = (
        system_residue_summary[system_residue_summary["system"] == system]
        .sort_values("weighted_mean_abs_delta_angstrom", ascending=False)
        .head(8)
    )
    clusters = sorted(set(edges["from_cluster"]).union(edges["to_cluster"]))
    inv = inventory[(inventory["system"] == system) & (inventory["cluster"].isin(clusters))]
    palette = ["red", "orange", "yelloworange", "tv_yellow", "salmon", "raspberry", "hotpink", "magenta"]
    system_color = "tv_blue" if system == "BM213" else "orange"

    lines = [
        f"# {system} major activation-core transition changes",
        "# Spheres highlight residues with largest transition-weighted changes.",
        "reinitialize",
        "bg_color white",
        "set cartoon_transparency, 0.25",
        "set sphere_quality, 2",
        "set label_size, 18",
    ]
    loaded_objects: list[str] = []
    for _, row in inv.sort_values("cluster").iterrows():
        if row["pdb_source"] == "missing" or not isinstance(row["pdb_path"], str):
            continue
        obj = f"{system}_C{int(row['cluster'])}"
        loaded_objects.append(obj)
        lines.extend(
            [
                f"load {row['pdb_path']}, {obj}",
                f"hide everything, {obj}",
                f"show cartoon, {obj}",
                f"color {system_color}, {obj}",
            ]
        )
    if loaded_objects:
        reference = loaded_objects[0]
        for obj in loaded_objects[1:]:
            lines.append(f"align {obj} and chain X and name CA, {reference} and chain X and name CA")

    for idx, (_, row) in enumerate(top.iterrows()):
        residue = row["residue"]
        resid = residue_to_canonical.get(residue)
        if resid is None:
            continue
        color = palette[min(idx, len(palette) - 1)]
        scale = 0.42 + min(float(row["weighted_mean_abs_delta_angstrom"]) / 1.4, 1.0) * 0.38
        sel = f"{system}_{residue}"
        lines.extend(
            [
                f"select {sel}, {' or '.join(loaded_objects)} and chain X and resi {resid} and name CA",
                f"show spheres, {sel}",
                f"set sphere_scale, {scale:.2f}, {sel}",
                f"color {color}, {sel}",
                f"label {sel}, \"{residue}\"",
            ]
        )
    lines.extend(["zoom", "orient"])
    out_pml.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_report(
    out_md: Path,
    dominant_edges: pd.DataFrame,
    residue_changes: pd.DataFrame,
    pair_changes: pd.DataFrame,
    system_residue_summary: pd.DataFrame,
    system_pair_summary: pd.DataFrame,
    displacement: pd.DataFrame,
    inventory: pd.DataFrame,
) -> None:
    lines: list[str] = []
    lines.append("# BM213/C5apep formal SOM activation-route transition changes")
    lines.append("")
    lines.append("## Method")
    lines.append("")
    lines.append(
        "- Used the 9-cluster formal transition network edge threshold shown in the figure: directed between-cluster edges with transitions >= 5."
    )
    lines.append(
        "- Quantified state-change residues from the 91 SOM activation-core distance features by comparing system-specific cluster mean distances before and after each transition."
    )
    lines.append(
        "- Residue magnitude is the mean absolute change across all activation-core distance features involving that residue; system-level magnitude is weighted by transition count."
    )
    lines.append(
        "- PDB support was computed by C-alpha aligning the representative structures and measuring per-residue C-alpha displacement for activation-core and AlloViz-route residues."
    )
    lines.append("")

    lines.append("## Dominant Transitions")
    lines.append("")
    for system in SYSTEMS:
        lines.append(f"### {system}")
        subset = dominant_edges[dominant_edges["system"] == system].sort_values(
            ["transitions", "from_cluster", "to_cluster"], ascending=[False, True, True]
        )
        for _, row in subset.iterrows():
            lines.append(
                f"- C{int(row['from_cluster'])} -> C{int(row['to_cluster'])}: "
                f"{int(row['transitions'])} transitions, P={row['transition_probability_from_cluster']:.3f}"
            )
        lines.append("")

    lines.append("## Main Residues By System")
    lines.append("")
    for system in SYSTEMS:
        lines.append(f"### {system}")
        subset = system_residue_summary[system_residue_summary["system"] == system].head(8)
        for _, row in subset.iterrows():
            lines.append(
                f"- {row['residue']}: weighted mean |delta|={row['weighted_mean_abs_delta_angstrom']:.2f} A; "
                f"max transition |delta|={row['max_mean_abs_delta_angstrom']:.2f} A; "
                f"top pairs: {row['top_pairs']}"
            )
        lines.append("")

    lines.append("## Main Residue Pairs By System")
    lines.append("")
    for system in SYSTEMS:
        lines.append(f"### {system}")
        subset = system_pair_summary[system_pair_summary["system"] == system].head(8)
        for _, row in subset.iterrows():
            direction = "increases" if row["weighted_signed_delta_angstrom"] > 0 else "decreases"
            lines.append(
                f"- {row['residue_pair']}: weighted |delta|={row['weighted_abs_delta_angstrom']:.2f} A; "
                f"weighted signed delta={row['weighted_signed_delta_angstrom']:.2f} A ({direction})"
            )
        lines.append("")

    lines.append("## Per-Transition Top Residues")
    lines.append("")
    for system in SYSTEMS:
        lines.append(f"### {system}")
        subset = dominant_edges[dominant_edges["system"] == system]
        for _, edge in subset.iterrows():
            trans_res = residue_changes[
                (residue_changes["system"] == system)
                & (residue_changes["from_cluster"] == edge["from_cluster"])
                & (residue_changes["to_cluster"] == edge["to_cluster"])
            ].head(5)
            residue_text = ", ".join(
                f"{r.residue} ({r.mean_abs_delta_angstrom:.2f} A)" for r in trans_res.itertuples()
            )
            top_pair_rows = pair_changes[
                (pair_changes["system"] == system)
                & (pair_changes["from_cluster"] == edge["from_cluster"])
                & (pair_changes["to_cluster"] == edge["to_cluster"])
            ].head(4)
            pair_text = ", ".join(
                f"{r.residue_pair} ({r.delta_angstrom:+.2f} A)" for r in top_pair_rows.itertuples()
            )
            lines.append(
                f"- C{int(edge['from_cluster'])} -> C{int(edge['to_cluster'])}: residues {residue_text}; pairs {pair_text}"
            )
        lines.append("")

    lines.append("## Representative PDB Sources")
    lines.append("")
    for system in SYSTEMS:
        lines.append(f"### {system}")
        subset = inventory[inventory["system"] == system].sort_values("cluster")
        for _, row in subset.iterrows():
            lines.append(f"- C{int(row['cluster'])}: {row['pdb_source']} - {row['pdb_name']}")
        lines.append("")

    if not displacement.empty:
        lines.append("## PDB C-alpha Displacement Support")
        lines.append("")
        for system in SYSTEMS:
            lines.append(f"### {system}")
            subset = displacement[displacement["system"] == system]
            for (from_c, to_c), group in subset.groupby(["from_cluster", "to_cluster"], sort=False):
                note = str(group["pdb_pair_note"].iloc[0]) if "pdb_pair_note" in group.columns else "ok"
                if note == "identical_ca_coordinates":
                    lines.append(
                        f"- C{int(from_c)} -> C{int(to_c)}: representative PDB C-alpha coordinates are identical; "
                        "PDB displacement is not informative for this edge."
                    )
                    continue
                top = group.sort_values("ca_displacement_angstrom", ascending=False).head(5)
                text = ", ".join(
                    f"{r.residue_label} ({r.ca_displacement_angstrom:.2f} A)" for r in top.itertuples()
                )
                lines.append(f"- C{int(from_c)} -> C{int(to_c)}: {text}")
            lines.append("")

    lines.append("## Notes")
    lines.append("")
    lines.append(
        "- These values describe structural/feature shifts associated with SOM state transitions; they are not by themselves causal proof of allosteric signal propagation."
    )
    lines.append(
        "- The PDB displacement table is a representative-structure check. The CSV feature tables are the primary quantitative basis because they average all frames assigned to each system-specific cluster."
    )
    out_md.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    change_dir = Path(__file__).resolve().parent
    formal_dir = change_dir.parent
    sommd_dir = formal_dir.parent.parent
    repo_dir = sommd_dir.parent.parent.parent
    input_dir = sommd_dir / "inputs" / "formal" / "activation_core"
    run_dir = sommd_dir / "runs" / "formal" / "activation_core"
    table_dir = sommd_dir / "tables"

    feature_matrix = pd.read_csv(input_dir / "feature_matrix.csv")
    frame_mapping = pd.read_csv(run_dir / "frame_mapping.csv")
    feature_meta = pd.read_csv(input_dir / "feature_metadata.csv")
    activation_mapping = pd.read_csv(input_dir / "activation_residue_mapping.csv")
    edges = pd.read_csv(table_dir / "formal_transition_edges.csv")

    route_summary_path = repo_dir / "article" / "analysis" / "alloviz" / "official_views" / "nglview" / "active" / "activation_routes_summary.csv"
    route_summary = pd.read_csv(route_summary_path) if route_summary_path.exists() else pd.DataFrame()

    if len(feature_matrix) != len(frame_mapping):
        raise ValueError("feature_matrix.csv and frame_mapping.csv have different row counts.")

    feature_cols = feature_matrix.columns.tolist()
    feature_info = feature_meta.set_index("feature_name").to_dict("index")
    frame_features = pd.concat(
        [frame_mapping[["global_frame", "system", "replicate", "cluster", "neuron"]], feature_matrix],
        axis=1,
    )
    frame_features = frame_features[frame_features["system"].isin(SYSTEMS)]

    cluster_means = frame_features.groupby(["system", "cluster"])[feature_cols].mean()
    cluster_frames = frame_features.groupby(["system", "cluster"]).size().rename("frames")

    dominant_edges = edges[
        (edges["system"].isin(SYSTEMS))
        & (edges["edge_type"] == "between")
        & (edges["transitions"] >= EDGE_THRESHOLD)
    ].copy()
    dominant_edges = dominant_edges.sort_values(
        ["system", "transitions", "from_cluster", "to_cluster"], ascending=[True, False, True, True]
    )

    pair_records: list[dict[str, object]] = []
    for edge in dominant_edges.itertuples(index=False):
        key_from = (edge.system, int(edge.from_cluster))
        key_to = (edge.system, int(edge.to_cluster))
        if key_from not in cluster_means.index or key_to not in cluster_means.index:
            continue
        from_vec = cluster_means.loc[key_from]
        to_vec = cluster_means.loc[key_to]
        for feature in feature_cols:
            info = feature_info[feature]
            delta = float(to_vec[feature] - from_vec[feature])
            pair_records.append(
                {
                    "system": edge.system,
                    "from_cluster": int(edge.from_cluster),
                    "to_cluster": int(edge.to_cluster),
                    "transitions": int(edge.transitions),
                    "transition_probability_from_cluster": float(edge.transition_probability_from_cluster),
                    "from_cluster_frames": int(cluster_frames.get(key_from, 0)),
                    "to_cluster_frames": int(cluster_frames.get(key_to, 0)),
                    "feature_name": feature,
                    "residue_a": info["residue_a"],
                    "residue_b": info["residue_b"],
                    "residue_pair": f"{info['residue_a']}-{info['residue_b']}",
                    "region_a": info["region_a"],
                    "region_b": info["region_b"],
                    "from_mean_distance_angstrom": float(from_vec[feature]),
                    "to_mean_distance_angstrom": float(to_vec[feature]),
                    "delta_angstrom": delta,
                    "abs_delta_angstrom": abs(delta),
                }
            )

    pair_changes = pd.DataFrame(pair_records)
    pair_changes["rank_abs_delta_in_transition"] = (
        pair_changes.groupby(["system", "from_cluster", "to_cluster"])["abs_delta_angstrom"]
        .rank(method="first", ascending=False)
        .astype(int)
    )
    pair_changes = pair_changes.sort_values(
        ["system", "from_cluster", "to_cluster", "rank_abs_delta_in_transition"]
    )

    residue_records: list[dict[str, object]] = []
    for (system, from_cluster, to_cluster), group in pair_changes.groupby(["system", "from_cluster", "to_cluster"]):
        for residue in CORE_RESIDUES:
            involved = group[(group["residue_a"] == residue) | (group["residue_b"] == residue)].copy()
            if involved.empty:
                continue
            top_pairs = involved.sort_values("abs_delta_angstrom", ascending=False).head(3)
            residue_records.append(
                {
                    "system": system,
                    "from_cluster": int(from_cluster),
                    "to_cluster": int(to_cluster),
                    "transitions": int(group["transitions"].iloc[0]),
                    "transition_probability_from_cluster": float(group["transition_probability_from_cluster"].iloc[0]),
                    "residue": residue,
                    "sum_abs_delta_angstrom": float(involved["abs_delta_angstrom"].sum()),
                    "mean_abs_delta_angstrom": float(involved["abs_delta_angstrom"].mean()),
                    "max_abs_delta_angstrom": float(involved["abs_delta_angstrom"].max()),
                    "signed_delta_mean_angstrom": float(involved["delta_angstrom"].mean()),
                    "top_pairs": "; ".join(
                        f"{r.residue_pair}:{r.delta_angstrom:+.2f}A" for r in top_pairs.itertuples()
                    ),
                }
            )

    residue_changes = pd.DataFrame(residue_records)
    residue_changes["rank_residue_in_transition"] = (
        residue_changes.groupby(["system", "from_cluster", "to_cluster"])["mean_abs_delta_angstrom"]
        .rank(method="first", ascending=False)
        .astype(int)
    )
    residue_changes = residue_changes.sort_values(
        ["system", "from_cluster", "to_cluster", "rank_residue_in_transition"]
    )

    system_residue_records: list[dict[str, object]] = []
    for (system, residue), group in residue_changes.groupby(["system", "residue"]):
        weights = group["transitions"].astype(float)
        weighted_mean = float(np.average(group["mean_abs_delta_angstrom"], weights=weights))
        top_pairs = (
            group.sort_values("mean_abs_delta_angstrom", ascending=False)["top_pairs"].head(2).tolist()
        )
        system_residue_records.append(
            {
                "system": system,
                "residue": residue,
                "weighted_mean_abs_delta_angstrom": weighted_mean,
                "max_mean_abs_delta_angstrom": float(group["mean_abs_delta_angstrom"].max()),
                "transitions_considered": int(group["transitions"].sum()),
                "top_pairs": " | ".join(top_pairs),
            }
        )
    system_residue_summary = pd.DataFrame(system_residue_records).sort_values(
        ["system", "weighted_mean_abs_delta_angstrom"], ascending=[True, False]
    )

    system_pair_records: list[dict[str, object]] = []
    for (system, feature), group in pair_changes.groupby(["system", "feature_name"]):
        weights = group["transitions"].astype(float)
        first = group.iloc[0]
        system_pair_records.append(
            {
                "system": system,
                "feature_name": feature,
                "residue_pair": first["residue_pair"],
                "residue_a": first["residue_a"],
                "residue_b": first["residue_b"],
                "weighted_abs_delta_angstrom": float(np.average(group["abs_delta_angstrom"], weights=weights)),
                "weighted_signed_delta_angstrom": float(np.average(group["delta_angstrom"], weights=weights)),
                "max_abs_delta_angstrom": float(group["abs_delta_angstrom"].max()),
            }
        )
    system_pair_summary = pd.DataFrame(system_pair_records).sort_values(
        ["system", "weighted_abs_delta_angstrom"], ascending=[True, False]
    )

    clusters_by_system: dict[str, set[int]] = defaultdict(set)
    for edge in dominant_edges.itertuples(index=False):
        clusters_by_system[edge.system].update([int(edge.from_cluster), int(edge.to_cluster)])

    inventory_records: list[dict[str, object]] = []
    representative_paths: dict[tuple[str, int], Path] = {}
    for system in SYSTEMS:
        for cluster in sorted(clusters_by_system[system]):
            path, source = find_representative_pdb(formal_dir, system, cluster)
            if path is not None:
                representative_paths[(system, cluster)] = path
            inventory_records.append(
                {
                    "system": system,
                    "cluster": cluster,
                    "pdb_source": source,
                    "pdb_name": path.name if path else "missing",
                    "pdb_path": str(path) if path else "",
                }
            )
    inventory = pd.DataFrame(inventory_records)

    core_canonical = (
        activation_mapping[activation_mapping["system"].isin(SYSTEMS)]
        .drop_duplicates("label")
        .set_index("label")["canonical_resid"]
        .astype(int)
        .to_dict()
    )
    route_residue_records: list[dict[str, object]] = []
    for system in SYSTEMS:
        route_labels = parse_route_residues(route_summary, system)
        for label in sorted(set(CORE_RESIDUES).union(route_labels), key=lambda x: label_to_nominal_number(x) or 0):
            canonical = core_canonical.get(label, label_to_canonical_resid(label))
            if canonical is None:
                continue
            route_residue_records.append(
                {
                    "system": system,
                    "residue_label": label,
                    "canonical_resid": int(canonical),
                    "is_som_activation_core": label in CORE_RESIDUES,
                    "is_alloviz_route_residue": label in route_labels,
                }
            )
    route_residues = pd.DataFrame(route_residue_records)

    displacement_records: list[dict[str, object]] = []
    pdb_pair_qc_records: list[dict[str, object]] = []
    pdb_cache: dict[Path, dict[int, tuple[str, np.ndarray]]] = {}
    for edge in dominant_edges.itertuples(index=False):
        from_path = representative_paths.get((edge.system, int(edge.from_cluster)))
        to_path = representative_paths.get((edge.system, int(edge.to_cluster)))
        if from_path is None or to_path is None:
            continue
        from_coords = pdb_cache.setdefault(from_path, read_ca_coords(from_path))
        to_coords = pdb_cache.setdefault(to_path, read_ca_coords(to_path))
        shared = sorted(set(from_coords) & set(to_coords))
        if len(shared) < 20:
            continue
        mobile = np.vstack([from_coords[resi][1] for resi in shared])
        target = np.vstack([to_coords[resi][1] for resi in shared])
        aligned_mobile = kabsch_align(mobile, target)
        aligned_by_resid = {resi: aligned_mobile[idx] for idx, resi in enumerate(shared)}
        aligned_distances = np.linalg.norm(target - aligned_mobile, axis=1)
        raw_distances = np.linalg.norm(target - mobile, axis=1)
        aligned_ca_rmsd = float(np.sqrt(np.mean(aligned_distances**2)))
        raw_max_ca_distance = float(np.max(raw_distances))
        pdb_pair_note = "identical_ca_coordinates" if raw_max_ca_distance < 1e-6 else "ok"
        pdb_pair_qc_records.append(
            {
                "system": edge.system,
                "from_cluster": int(edge.from_cluster),
                "to_cluster": int(edge.to_cluster),
                "transitions": int(edge.transitions),
                "shared_ca_residues": len(shared),
                "aligned_ca_rmsd_angstrom": aligned_ca_rmsd,
                "raw_max_ca_distance_angstrom": raw_max_ca_distance,
                "pdb_pair_note": pdb_pair_note,
                "from_pdb": str(from_path),
                "to_pdb": str(to_path),
            }
        )
        route_subset = route_residues[route_residues["system"] == edge.system]
        for route_row in route_subset.itertuples(index=False):
            resid = int(route_row.canonical_resid)
            if resid not in from_coords or resid not in to_coords:
                continue
            displacement = float(np.linalg.norm(to_coords[resid][1] - aligned_by_resid[resid]))
            displacement_records.append(
                {
                    "system": edge.system,
                    "from_cluster": int(edge.from_cluster),
                    "to_cluster": int(edge.to_cluster),
                    "transitions": int(edge.transitions),
                    "residue_label": route_row.residue_label,
                    "canonical_resid": resid,
                    "from_resname": from_coords[resid][0],
                    "to_resname": to_coords[resid][0],
                    "ca_displacement_angstrom": displacement,
                    "aligned_ca_rmsd_all_shared_angstrom": aligned_ca_rmsd,
                    "pdb_pair_note": pdb_pair_note,
                    "is_som_activation_core": bool(route_row.is_som_activation_core),
                    "is_alloviz_route_residue": bool(route_row.is_alloviz_route_residue),
                    "from_pdb": str(from_path),
                    "to_pdb": str(to_path),
                }
            )
    displacement = pd.DataFrame(displacement_records)
    if not displacement.empty:
        displacement["rank_displacement_in_transition"] = (
            displacement.groupby(["system", "from_cluster", "to_cluster"])["ca_displacement_angstrom"]
            .rank(method="first", ascending=False)
            .astype(int)
        )
        displacement = displacement.sort_values(
            ["system", "from_cluster", "to_cluster", "rank_displacement_in_transition"]
        )
    pdb_pair_qc = pd.DataFrame(pdb_pair_qc_records)

    pair_changes.to_csv(change_dir / "transition_activation_pair_changes.csv", index=False)
    residue_changes.to_csv(change_dir / "transition_activation_residue_changes.csv", index=False)
    system_residue_summary.to_csv(change_dir / "system_activation_residue_change_summary.csv", index=False)
    system_pair_summary.to_csv(change_dir / "system_activation_pair_change_summary.csv", index=False)
    dominant_edges.to_csv(change_dir / "dominant_transition_edges_used.csv", index=False)
    inventory.to_csv(change_dir / "pdb_representative_inventory.csv", index=False)
    route_residues.to_csv(change_dir / "route_residues_for_pdb_displacement.csv", index=False)
    displacement.to_csv(change_dir / "pdb_ca_displacement_by_transition.csv", index=False)
    pdb_pair_qc.to_csv(change_dir / "pdb_pair_qc.csv", index=False)

    make_heatmap(residue_changes, change_dir / "activation_transition_residue_heatmap.png")
    make_heatmap(
        residue_changes,
        change_dir / "activation_transition_residue_heatmap_no_panel_labels.png",
        show_panel_labels=False,
    )
    make_summary_barplot(system_residue_summary, change_dir / "activation_transition_residue_summary.png")

    residue_to_canonical = {label: int(resid) for label, resid in core_canonical.items()}
    for system in SYSTEMS:
        write_pml(
            system=system,
            edges=dominant_edges[dominant_edges["system"] == system],
            inventory=inventory,
            system_residue_summary=system_residue_summary,
            residue_to_canonical=residue_to_canonical,
            out_pml=change_dir / f"{system}_activation_transition_changes.pml",
        )

    write_report(
        change_dir / "transition_activation_change_summary.md",
        dominant_edges,
        residue_changes,
        pair_changes,
        system_residue_summary,
        system_pair_summary,
        displacement,
        inventory,
    )

    print(f"Wrote transition-change outputs to {change_dir}")
    print("Dominant edges:")
    print(dominant_edges[["system", "from_cluster", "to_cluster", "transitions"]].to_string(index=False))
    print("\nTop residues:")
    print(system_residue_summary.groupby("system").head(6).to_string(index=False))


if __name__ == "__main__":
    main()
