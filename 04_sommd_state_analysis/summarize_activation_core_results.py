#!/usr/bin/env python3
from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
import MDAnalysis as mda
import numpy as np
import pandas as pd


ROOT = Path("/mnt/e/work/modeling/c5ar/md")
SOMMD_ROOT = ROOT / "article/analysis/sommd"
SYSTEM_ORDER = ("apo", "BM213", "C5apep")
SYSTEM_COLORS = {
    "apo": "#7a7a7a",
    "BM213": "#1f78b4",
    "C5apep": "#e68613",
    "all": "#4d4d4d",
}
COMPARISONS = (
    ("BM213", "apo"),
    ("C5apep", "apo"),
    ("C5apep", "BM213"),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Summarize C5aR1 activation-core SOM outputs.")
    parser.add_argument("--profile", choices=("smoke", "pilot", "formal"), required=True)
    parser.add_argument(
        "--run-dir",
        default=None,
        help="SOM run directory. Defaults to article/analysis/sommd/runs/<profile>/activation_core.",
    )
    parser.add_argument(
        "--input-dir",
        default=None,
        help="Feature input directory. Defaults to article/analysis/sommd/inputs/<profile>/activation_core.",
    )
    return parser.parse_args()


def resolve_paths(args: argparse.Namespace) -> tuple[Path, Path]:
    run_dir = Path(args.run_dir) if args.run_dir else SOMMD_ROOT / "runs" / args.profile / "activation_core"
    input_dir = Path(args.input_dir) if args.input_dir else SOMMD_ROOT / "inputs" / args.profile / "activation_core"
    if not run_dir.is_absolute():
        run_dir = ROOT / run_dir
    if not input_dir.is_absolute():
        input_dir = ROOT / input_dir
    return run_dir, input_dir


def ensure_dirs(profile: str) -> dict[str, Path]:
    dirs = {
        "tables": SOMMD_ROOT / "tables",
        "figures": SOMMD_ROOT / "figures",
        "pymol": SOMMD_ROOT / "pymol" / profile,
    }
    for path in dirs.values():
        path.mkdir(parents=True, exist_ok=True)
    return dirs


def read_inputs(run_dir: Path, input_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    frame_map = pd.read_csv(run_dir / "frame_mapping.csv")
    feature_matrix = pd.read_csv(input_dir / "feature_matrix.csv")
    feature_metadata = pd.read_csv(input_dir / "feature_metadata.csv")
    residue_mapping = pd.read_csv(input_dir / "activation_residue_mapping.csv")
    if len(frame_map) != len(feature_matrix):
        raise ValueError("frame_mapping and feature_matrix row counts differ")
    return frame_map, feature_matrix, feature_metadata, residue_mapping


def population_tables(frame_map: pd.DataFrame, out_tables: Path, profile: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    clusters = sorted(frame_map["cluster"].astype(int).unique())
    system_counts = (
        frame_map.groupby(["system", "cluster"], observed=True)
        .size()
        .rename("frames")
        .reset_index()
    )
    all_rows = pd.MultiIndex.from_product([SYSTEM_ORDER, clusters], names=["system", "cluster"]).to_frame(index=False)
    system_counts = all_rows.merge(system_counts, on=["system", "cluster"], how="left").fillna({"frames": 0})
    system_counts["frames"] = system_counts["frames"].astype(int)
    total_by_system = system_counts.groupby("system")["frames"].transform("sum")
    system_counts["fraction_within_system"] = system_counts["frames"] / total_by_system
    system_counts.to_csv(out_tables / f"{profile}_cluster_population_by_system.csv", index=False)

    rep_counts = (
        frame_map.groupby(["system", "replicate", "cluster"], observed=True)
        .size()
        .rename("frames")
        .reset_index()
    )
    reps = frame_map[["system", "replicate"]].drop_duplicates()
    rep_grid = reps.merge(pd.DataFrame({"cluster": clusters}), how="cross")
    rep_counts = rep_grid.merge(rep_counts, on=["system", "replicate", "cluster"], how="left").fillna({"frames": 0})
    rep_counts["frames"] = rep_counts["frames"].astype(int)
    total_by_rep = rep_counts.groupby(["system", "replicate"])["frames"].transform("sum")
    rep_counts["fraction_within_replicate"] = rep_counts["frames"] / total_by_rep
    rep_counts.to_csv(out_tables / f"{profile}_cluster_population_by_replicate.csv", index=False)
    return system_counts, rep_counts


def enrichment_table(system_counts: pd.DataFrame, out_tables: Path, profile: str) -> pd.DataFrame:
    pivot = system_counts.pivot(index="cluster", columns="system", values="fraction_within_system").fillna(0.0)
    rows: list[dict[str, object]] = []
    eps = 1e-6
    for cluster, values in pivot.iterrows():
        for left, right in COMPARISONS:
            left_fraction = float(values.get(left, 0.0))
            right_fraction = float(values.get(right, 0.0))
            rows.append(
                {
                    "cluster": int(cluster),
                    "comparison": f"{left}_minus_{right}",
                    "left_system": left,
                    "right_system": right,
                    "left_fraction": left_fraction,
                    "right_fraction": right_fraction,
                    "fraction_delta": left_fraction - right_fraction,
                    "log2_fraction_ratio": float(np.log2((left_fraction + eps) / (right_fraction + eps))),
                }
            )
    out = pd.DataFrame(rows)
    out.to_csv(out_tables / f"{profile}_cluster_enrichment_vs_apo.csv", index=False)
    return out


def feature_by_cluster(
    frame_map: pd.DataFrame,
    feature_matrix: pd.DataFrame,
    feature_metadata: pd.DataFrame,
    out_tables: Path,
    profile: str,
) -> pd.DataFrame:
    features = feature_metadata["feature_name"].tolist()
    data = pd.concat([frame_map[["cluster"]].reset_index(drop=True), feature_matrix[features].reset_index(drop=True)], axis=1)
    means = data.groupby("cluster")[features].mean()
    overall = feature_matrix[features].mean()
    sd = feature_matrix[features].std(ddof=1).replace(0, 1.0)
    rows: list[dict[str, object]] = []
    meta = feature_metadata.set_index("feature_name")
    for cluster, values in means.iterrows():
        for feature in features:
            m = meta.loc[feature]
            rows.append(
                {
                    "cluster": int(cluster),
                    "feature_name": feature,
                    "mean_distance_angstrom": float(values[feature]),
                    "overall_mean_distance_angstrom": float(overall[feature]),
                    "delta_from_overall_angstrom": float(values[feature] - overall[feature]),
                    "z_delta_from_overall": float((values[feature] - overall[feature]) / sd[feature]),
                    "residue_a": m["residue_a"],
                    "residue_b": m["residue_b"],
                    "region_a": m["region_a"],
                    "region_b": m["region_b"],
                }
            )
    out = pd.DataFrame(rows)
    out.to_csv(out_tables / f"{profile}_activation_feature_by_cluster.csv", index=False)
    return out


def plot_population(system_counts: pd.DataFrame, out_figures: Path, profile: str) -> Path:
    pivot = system_counts.pivot(index="system", columns="cluster", values="fraction_within_system")
    pivot = pivot.reindex(SYSTEM_ORDER)
    fig, ax = plt.subplots(figsize=(9, 3.2))
    im = ax.imshow(pivot.values, aspect="auto", cmap="viridis", vmin=0)
    ax.set_xticks(range(len(pivot.columns)))
    ax.set_xticklabels([str(c) for c in pivot.columns])
    ax.set_yticks(range(len(pivot.index)))
    ax.set_yticklabels(pivot.index)
    ax.set_xlabel("SOM cluster")
    ax.set_title(f"{profile} activation-core cluster occupancy")
    fig.colorbar(im, ax=ax, label="fraction within system")
    fig.tight_layout()
    out = out_figures / f"{profile}_cluster_population_heatmap.png"
    fig.savefig(out, dpi=300)
    plt.close(fig)
    return out


def plot_feature_heatmap(feature_table: pd.DataFrame, feature_metadata: pd.DataFrame, out_figures: Path, profile: str) -> Path:
    feature_order = feature_metadata["feature_name"].tolist()
    pivot = feature_table.pivot(index="feature_name", columns="cluster", values="z_delta_from_overall")
    pivot = pivot.reindex(feature_order)
    labels = [
        f"{row.residue_a}-{row.residue_b}"
        for row in feature_metadata.itertuples(index=False)
    ]
    height = max(12.0, 0.16 * len(feature_order))
    fig, ax = plt.subplots(figsize=(8, height))
    vmax = max(1.0, float(np.nanmax(np.abs(pivot.values))))
    im = ax.imshow(pivot.values, aspect="auto", cmap="coolwarm", vmin=-vmax, vmax=vmax)
    ax.set_xticks(range(len(pivot.columns)))
    ax.set_xticklabels([str(c) for c in pivot.columns])
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=5)
    ax.set_xlabel("SOM cluster")
    ax.set_title(f"{profile} activation-core feature shifts")
    fig.colorbar(im, ax=ax, label="cluster mean z-delta")
    fig.tight_layout()
    out = out_figures / f"{profile}_activation_feature_heatmap.png"
    fig.savefig(out, dpi=300)
    plt.close(fig)
    return out


def plot_enrichment(enrichment: pd.DataFrame, out_figures: Path, profile: str) -> Path:
    comparisons = enrichment["comparison"].drop_duplicates().tolist()
    clusters = sorted(enrichment["cluster"].unique())
    fig, axes = plt.subplots(len(comparisons), 1, figsize=(8, 2.6 * len(comparisons)), sharex=True)
    if len(comparisons) == 1:
        axes = [axes]
    for ax, comparison in zip(axes, comparisons):
        sub = enrichment[enrichment["comparison"] == comparison].set_index("cluster").reindex(clusters)
        ax.axhline(0, color="#333333", lw=0.8)
        colors = ["#1f78b4" if x >= 0 else "#b15928" for x in sub["fraction_delta"]]
        ax.bar(clusters, sub["fraction_delta"], color=colors)
        ax.set_ylabel("fraction delta")
        ax.set_title(comparison)
    axes[-1].set_xlabel("SOM cluster")
    fig.tight_layout()
    out = out_figures / f"{profile}_cluster_enrichment.png"
    fig.savefig(out, dpi=300)
    plt.close(fig)
    return out


def transition_tables(frame_map: pd.DataFrame, out_tables: Path, profile: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    clusters = sorted(frame_map["cluster"].astype(int).unique())
    cluster_to_idx = {cluster: idx for idx, cluster in enumerate(clusters)}
    systems = ("all",) + SYSTEM_ORDER
    matrix_rows: list[dict[str, object]] = []
    edge_rows: list[dict[str, object]] = []
    node_rows: list[dict[str, object]] = []

    def count_transitions(sub: pd.DataFrame) -> np.ndarray:
        matrix = np.zeros((len(clusters), len(clusters)), dtype=int)
        for _, seg in sub.groupby(["system", "replicate"], observed=True, sort=False):
            seg = seg.sort_values(["local_frame_index", "som_frame"])
            values = seg["cluster"].astype(int).to_numpy()
            if len(values) < 2:
                continue
            for left, right in zip(values[:-1], values[1:]):
                matrix[cluster_to_idx[int(left)], cluster_to_idx[int(right)]] += 1
        return matrix

    for system in systems:
        sub = frame_map if system == "all" else frame_map[frame_map["system"] == system]
        matrix = count_transitions(sub)
        population = sub["cluster"].astype(int).value_counts().reindex(clusters, fill_value=0)
        total_frames = int(population.sum())
        total_transitions = int(matrix.sum())
        for cluster in clusters:
            frames = int(population.loc[cluster])
            node_rows.append(
                {
                    "profile": profile,
                    "system": system,
                    "cluster": cluster,
                    "frames": frames,
                    "fraction_within_system": frames / total_frames if total_frames else 0.0,
                }
            )
        for i, from_cluster in enumerate(clusters):
            from_total = int(matrix[i, :].sum())
            for j, to_cluster in enumerate(clusters):
                transitions = int(matrix[i, j])
                row = {
                    "profile": profile,
                    "system": system,
                    "from_cluster": from_cluster,
                    "to_cluster": to_cluster,
                    "transitions": transitions,
                    "from_total_transitions": from_total,
                    "transition_probability_from_cluster": transitions / from_total if from_total else 0.0,
                    "fraction_of_system_transitions": transitions / total_transitions if total_transitions else 0.0,
                    "edge_type": "self" if from_cluster == to_cluster else "between",
                }
                matrix_rows.append(row)
                if transitions > 0:
                    edge_rows.append(row)

    matrix_table = pd.DataFrame(matrix_rows)
    edge_table = pd.DataFrame(edge_rows)
    node_table = pd.DataFrame(node_rows)
    matrix_table.to_csv(out_tables / f"{profile}_transition_matrix_long.csv", index=False)
    edge_table.to_csv(out_tables / f"{profile}_transition_edges.csv", index=False)
    node_table.to_csv(out_tables / f"{profile}_transition_nodes.csv", index=False)
    return matrix_table, edge_table, node_table


def cluster_positions(clusters: list[int]) -> dict[int, tuple[float, float]]:
    angles = np.linspace(np.pi / 2.0, np.pi / 2.0 - 2.0 * np.pi, len(clusters), endpoint=False)
    return {
        cluster: (float(np.cos(angle)), float(np.sin(angle)))
        for cluster, angle in zip(clusters, angles)
    }


def draw_transition_network(
    ax: plt.Axes,
    clusters: list[int],
    edges: pd.DataFrame,
    nodes: pd.DataFrame,
    system: str,
    title: str,
    min_count: int,
) -> None:
    positions = cluster_positions(clusters)
    sub_edges = edges[
        (edges["system"] == system)
        & (edges["edge_type"] == "between")
        & (edges["transitions"] >= min_count)
    ].copy()
    sub_nodes = nodes[nodes["system"] == system].set_index("cluster")
    max_edge = int(sub_edges["transitions"].max()) if not sub_edges.empty else 1

    for row in sub_edges.itertuples(index=False):
        x1, y1 = positions[int(row.from_cluster)]
        x2, y2 = positions[int(row.to_cluster)]
        ratio = float(row.transitions) / max_edge
        rad = 0.18 if int(row.from_cluster) < int(row.to_cluster) else -0.18
        arrow = FancyArrowPatch(
            (x1, y1),
            (x2, y2),
            arrowstyle="-|>",
            mutation_scale=8.0 + 9.0 * ratio,
            linewidth=0.7 + 4.2 * (ratio ** 0.7),
            color="#404040",
            alpha=0.22 + 0.62 * ratio,
            connectionstyle=f"arc3,rad={rad}",
            shrinkA=20,
            shrinkB=20,
            zorder=1,
        )
        ax.add_patch(arrow)

    label_edges = sub_edges.sort_values("transitions", ascending=False).head(5)
    for row in label_edges.itertuples(index=False):
        x1, y1 = positions[int(row.from_cluster)]
        x2, y2 = positions[int(row.to_cluster)]
        mx = (x1 + x2) / 2.0
        my = (y1 + y2) / 2.0
        ax.text(
            mx,
            my,
            str(int(row.transitions)),
            ha="center",
            va="center",
            fontsize=7,
            color="#222222",
            bbox={"boxstyle": "round,pad=0.18", "facecolor": "white", "edgecolor": "none", "alpha": 0.75},
            zorder=4,
        )

    for cluster in clusters:
        x, y = positions[cluster]
        fraction = float(sub_nodes.loc[cluster, "fraction_within_system"]) if cluster in sub_nodes.index else 0.0
        size = 260.0 + 2600.0 * np.sqrt(max(fraction, 0.0))
        alpha = 0.95 if fraction > 0 else 0.32
        ax.scatter(
            [x],
            [y],
            s=size,
            color=SYSTEM_COLORS.get(system, "#4d4d4d"),
            alpha=alpha,
            edgecolor="white",
            linewidth=1.6,
            zorder=3,
        )
        ax.text(x, y, str(cluster), ha="center", va="center", fontsize=9, color="white", weight="bold", zorder=5)

    ax.set_title(title, fontsize=11)
    ax.text(
        0.0,
        -1.35,
        f"edges >= {min_count}; node size = occupancy",
        ha="center",
        va="center",
        fontsize=8,
        color="#555555",
    )
    ax.set_xlim(-1.45, 1.45)
    ax.set_ylim(-1.45, 1.45)
    ax.set_aspect("equal")
    ax.axis("off")


def plot_transition_networks(
    edges: pd.DataFrame,
    nodes: pd.DataFrame,
    out_figures: Path,
    profile: str,
    min_count: int = 5,
) -> list[Path]:
    clusters = sorted(nodes["cluster"].astype(int).unique())
    out_paths: list[Path] = []

    fig, ax = plt.subplots(figsize=(6.2, 5.7))
    draw_transition_network(ax, clusters, edges, nodes, "all", f"{profile} activation-core transition network", min_count)
    fig.tight_layout()
    out = out_figures / f"{profile}_transition_network_overall.png"
    fig.savefig(out, dpi=300)
    plt.close(fig)
    out_paths.append(out)

    fig, axes = plt.subplots(1, len(SYSTEM_ORDER), figsize=(15.0, 5.2))
    for ax, system in zip(axes, SYSTEM_ORDER):
        draw_transition_network(ax, clusters, edges, nodes, system, system, min_count)
    fig.suptitle(f"{profile} activation-core transition networks by system", fontsize=13)
    fig.tight_layout()
    out = out_figures / f"{profile}_transition_network_by_system.png"
    fig.savefig(out, dpi=300)
    plt.close(fig)
    out_paths.append(out)
    return out_paths


def export_representatives(
    run_dir: Path,
    frame_map: pd.DataFrame,
    residue_mapping: pd.DataFrame,
    out_tables: Path,
    out_pymol: Path,
    profile: str,
) -> pd.DataFrame:
    reps = pd.read_csv(run_dir / "representative_frames.csv")
    rows: list[dict[str, object]] = []
    for rep in reps.itertuples(index=False):
        som_frame = int(rep.representative_frame)
        hit = frame_map[frame_map["som_frame"].astype(int) == som_frame]
        if len(hit) != 1:
            raise ValueError(f"Representative frame {som_frame}: expected one frame metadata row, found {len(hit)}")
        meta = hit.iloc[0]
        pdb = ROOT / str(meta["pdb"])
        traj = ROOT / str(meta["downsampled_trajectory"])
        local_frame_index = int(meta["local_frame_index"])
        universe = mda.Universe(str(pdb), str(traj))
        universe.trajectory[local_frame_index - 1]
        object_name = f"cluster_{int(rep.cluster):02d}_{meta['system']}_{meta['replicate']}"
        out_pdb = out_pymol / f"{object_name}.pdb"
        universe.atoms.write(str(out_pdb))
        rows.append(
            {
                "cluster": int(rep.cluster),
                "representative_som_frame": som_frame,
                "representative_neuron": int(rep.representative_neuron),
                "distance_to_neuron": float(rep.distance_to_neuron),
                "system": meta["system"],
                "replicate": meta["replicate"],
                "local_frame_index": local_frame_index,
                "time_ns": float(meta["time_ns"]),
                "pdb": meta["pdb"],
                "downsampled_trajectory": meta["downsampled_trajectory"],
                "representative_pdb": out_pdb.relative_to(ROOT).as_posix(),
                "object_name": object_name,
            }
        )
    out = pd.DataFrame(rows)
    out.to_csv(out_tables / f"{profile}_representative_frames_full.csv", index=False)
    write_pymol_script(out, residue_mapping, out_pymol, profile)
    return out


def write_pymol_script(representatives: pd.DataFrame, residue_mapping: pd.DataFrame, out_pymol: Path, profile: str) -> Path:
    colors = {"apo": "gray70", "BM213": "tv_blue", "C5apep": "orange"}
    lines = [
        "reinitialize",
        "bg_color white",
        "set cartoon_transparency, 0.15",
        "set sphere_scale, 0.45",
        "set label_size, 14",
    ]
    for row in representatives.itertuples(index=False):
        pdb_path = ROOT / row.representative_pdb
        obj = row.object_name
        system = row.system
        raw = residue_mapping[residue_mapping["system"] == system]["raw_resid"].astype(int).tolist()
        resid_expr = "+".join(str(x) for x in raw)
        lines.extend(
            [
                f"load {pdb_path.as_posix()}, {obj}",
                f"hide everything, {obj}",
                f"show cartoon, {obj}",
                f"color {colors.get(system, 'gray50')}, {obj}",
                f"select activation_core_{obj}, {obj} and resi {resid_expr} and name CA",
                f"show spheres, activation_core_{obj}",
                f"color yellow, activation_core_{obj}",
                f"label activation_core_{obj}, resn + resi",
            ]
        )
    if not representatives.empty:
        first = representatives.iloc[0]["object_name"]
        for obj in representatives["object_name"].iloc[1:]:
            lines.append(f"align {obj} and name CA, {first} and name CA")
    lines.extend(["orient", "zoom", "set ray_trace_mode, 1"])
    out = out_pymol / f"{profile}_load_cluster_representatives.pml"
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return out


def write_report(
    profile: str,
    run_dir: Path,
    input_dir: Path,
    tables: dict[str, Path],
    figures: list[Path],
    representatives: pd.DataFrame,
    out_tables: Path,
) -> Path:
    report = out_tables.parent / f"report_{profile}.md"
    lines = [
        f"# C5aR1 SOMMD Activation-Core Report ({profile})",
        "",
        "## Inputs",
        "",
        f"- Feature input: `{input_dir.relative_to(ROOT).as_posix()}`",
        f"- SOM run: `{run_dir.relative_to(ROOT).as_posix()}`",
        "- Feature definition: 14 activation-core residues, 91 CA-CA distance features.",
        "- Main interpretation: apo = basal reference; BM213 = G protein-compatible partial agonist state; C5apep = broader full agonist ensemble.",
        "",
        "## Tables",
        "",
    ]
    for label, path in tables.items():
        lines.append(f"- {label}: `{path.relative_to(ROOT).as_posix()}`")
    lines.extend(["", "## Figures", ""])
    for path in figures:
        lines.append(f"- `{path.relative_to(ROOT).as_posix()}`")
    lines.extend(["", "## Representative Structures", ""])
    for row in representatives.itertuples(index=False):
        lines.append(
            f"- cluster {row.cluster}: {row.system} {row.replicate}, frame {row.local_frame_index}, "
            f"`{row.representative_pdb}`"
        )
    lines.extend(
        [
            "",
            "## Interpretation Guide",
            "",
            "- apo-enriched clusters are basal or inactive-like candidates.",
            "- BM213-enriched clusters are G protein-compatible activation candidates.",
            "- C5apep clusters overlapping BM213 support G protein-compatible states.",
            "- C5apep-specific enriched clusters are candidates for broader activation or beta-arrestin-compatible states.",
        ]
    )
    report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    if profile == "formal":
        shutil.copyfile(report, out_tables.parent / "report.md")
    return report


def main() -> int:
    args = parse_args()
    run_dir, input_dir = resolve_paths(args)
    dirs = ensure_dirs(args.profile)
    frame_map, feature_matrix, feature_metadata, residue_mapping = read_inputs(run_dir, input_dir)

    system_counts, rep_counts = population_tables(frame_map, dirs["tables"], args.profile)
    enrichment = enrichment_table(system_counts, dirs["tables"], args.profile)
    feature_table = feature_by_cluster(frame_map, feature_matrix, feature_metadata, dirs["tables"], args.profile)
    representatives = export_representatives(
        run_dir,
        frame_map,
        residue_mapping,
        dirs["tables"],
        dirs["pymol"],
        args.profile,
    )
    figures = [
        plot_population(system_counts, dirs["figures"], args.profile),
        plot_feature_heatmap(feature_table, feature_metadata, dirs["figures"], args.profile),
        plot_enrichment(enrichment, dirs["figures"], args.profile),
    ]
    table_paths = {
        "cluster_population_by_system": dirs["tables"] / f"{args.profile}_cluster_population_by_system.csv",
        "cluster_population_by_replicate": dirs["tables"] / f"{args.profile}_cluster_population_by_replicate.csv",
        "cluster_enrichment": dirs["tables"] / f"{args.profile}_cluster_enrichment_vs_apo.csv",
        "activation_feature_by_cluster": dirs["tables"] / f"{args.profile}_activation_feature_by_cluster.csv",
        "representative_frames_full": dirs["tables"] / f"{args.profile}_representative_frames_full.csv",
    }
    if args.profile == "formal":
        transition_matrix, transition_edges, transition_nodes = transition_tables(frame_map, dirs["tables"], args.profile)
        figures.extend(plot_transition_networks(transition_edges, transition_nodes, dirs["figures"], args.profile))
        table_paths.update(
            {
                "transition_matrix_long": dirs["tables"] / f"{args.profile}_transition_matrix_long.csv",
                "transition_edges": dirs["tables"] / f"{args.profile}_transition_edges.csv",
                "transition_nodes": dirs["tables"] / f"{args.profile}_transition_nodes.csv",
            }
        )
    report = write_report(args.profile, run_dir, input_dir, table_paths, figures, representatives, dirs["tables"])
    print(report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
