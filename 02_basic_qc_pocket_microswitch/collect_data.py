from __future__ import annotations

import ast
import re
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.vq import kmeans2, vq

from figure_utils import (
    DATA_DIR,
    FIG_ROOT,
    LOG_DIR,
    METRIC_FILES,
    ROOT,
    SYSTEMS,
    TABLE_DIR,
    add_time_ns,
    ensure_dirs,
    metric_path,
    note_warning,
    parse_production_log,
    read_dat,
    summarize_series,
    write_placeholder_readme,
)


HBOND_METRICS = ["Y258_hbond", "S171_hbond"]
MICROSWITCH_METRICS = [
    "M120_distance",
    "W255_chi2",
    "R350_chi2",
    "TM62_IC_distance",
    "TM72_IC_distance",
    "DRY_PIF_distance",
]
RMSD_METRICS = ["RMSD"]


def collect_metric(metric: str, warnings: list[str]) -> pd.DataFrame:
    frames = []
    for info in SYSTEMS:
        system = info["system"]
        rep = info["replicate"]
        key = (system, rep, metric)
        if key not in METRIC_FILES:
            warnings.append(f"Missing metric mapping: {key}")
            continue
        path = metric_path(system, rep, metric)
        if not path.exists():
            warnings.append(f"Missing metric file: {path.relative_to(ROOT)}")
            continue
        df = read_dat(path)
        if df.empty:
            warnings.append(f"Empty or unparsable metric file: {path.relative_to(ROOT)}")
            continue
        df = add_time_ns(df)
        df["system"] = system
        df["replicate"] = rep
        df["metric"] = metric
        df["source_path"] = str(path.relative_to(ROOT))
        frames.append(df[["system", "replicate", "metric", "frame", "time_ns", "value", "source_path"]])
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def collect_inventory(all_series: list[pd.DataFrame]) -> pd.DataFrame:
    rows = []
    for df in all_series:
        if df.empty:
            continue
        for (system, rep, metric, src), sub in df.groupby(["system", "replicate", "metric", "source_path"]):
            rows.append(
                {
                    "system": system,
                    "replicate": rep,
                    "metric": metric,
                    "source_path": src,
                    "frames": int(len(sub)),
                    "time_ns_min": round(float(sub["time_ns"].min()), 3),
                    "time_ns_max": round(float(sub["time_ns"].max()), 3),
                    "value_mean": round(float(sub["value"].mean()), 5),
                    "value_min": round(float(sub["value"].min()), 5),
                    "value_max": round(float(sub["value"].max()), 5),
                }
            )
    return pd.DataFrame(rows)


def collect_qc() -> pd.DataFrame:
    patterns = {
        ("BM213", "rep1"): ["production*.out"],
        ("BM213", "rep2"): ["production*.out"],
        ("BM213", "rep3"): ["production*.out"],
        ("C5apep", "rep1"): ["production*.out"],
        ("C5apep", "rep2"): ["production*.out"],
        ("C5apep", "rep3"): ["production*.out"],
    }
    rows = []
    for info in SYSTEMS:
        system = info["system"]
        rep = info["replicate"]
        files = []
        for pattern in patterns[(system, rep)]:
            files.extend(sorted(info["dir"].glob(pattern)))
        for path in files:
            rec = parse_production_log(path)
            rec["system"] = system
            rec["replicate"] = rep
            rows.append(rec)
    cols = ["system", "replicate", "log_path", "exists", "nstlim", "dt_ps", "temp0_k", "ntwx", "last_nstep", "last_time_ps", "mean_temp_k", "mean_density", "error_or_warning_count"]
    return pd.DataFrame(rows, columns=cols)


def make_fes_points(microswitch: pd.DataFrame, warnings: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    for system in ["BM213", "C5apep"]:
        for rep in ["rep1", "rep2", "rep3"]:
            x = microswitch[(microswitch["system"] == system) & (microswitch["replicate"] == rep) & (microswitch["metric"] == "M120_distance")]
            y = microswitch[(microswitch["system"] == system) & (microswitch["replicate"] == rep) & (microswitch["metric"] == "DRY_PIF_distance")]
            n = min(len(x), len(y))
            if n == 0:
                warnings.append(f"Cannot build FES points for {system} {rep}; missing M120 or DRY/PIF data")
                continue
            if len(x) != len(y):
                warnings.append(f"FES frame mismatch for {system} {rep}: M120={len(x)}, DRY/PIF={len(y)}; truncated to {n}")
            x = x.sort_values("frame").iloc[:n].reset_index(drop=True)
            y = y.sort_values("frame").iloc[:n].reset_index(drop=True)
            rows.append(
                pd.DataFrame(
                    {
                        "system": system,
                        "replicate": rep,
                        "frame": x["frame"],
                        "time_ns": x["time_ns"],
                        "M120_distance": x["value"],
                        "DRY_PIF_distance": y["value"],
                        "source_M120": x["source_path"],
                        "source_DRY_PIF": y["source_path"],
                    }
                )
            )
    fes = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    if fes.empty:
        return fes, pd.DataFrame()

    xy = fes[["M120_distance", "DRY_PIF_distance"]].to_numpy(float)
    finite = np.isfinite(xy).all(axis=1)
    xy = xy[finite]
    mean = xy.mean(axis=0)
    std = xy.std(axis=0)
    std[std == 0] = 1
    z = (xy - mean) / std
    sample_idx = np.linspace(0, len(z) - 1, min(12000, len(z))).astype(int)
    centroids_z, _ = kmeans2(z[sample_idx], 3, minit="points", iter=60, seed=7)
    labels, _ = vq(z, centroids_z)
    centroids = centroids_z * std + mean
    order = np.argsort(centroids[:, 0])
    label_map = {old: f"Basin {i + 1}" for i, old in enumerate(order)}
    centroids_map = {label_map[i]: centroids[i] for i in range(len(centroids))}
    fes = fes.loc[finite].copy()
    fes["basin"] = [label_map[int(x)] for x in labels]

    occ = (
        fes.groupby(["system", "replicate", "basin"], as_index=False)
        .size()
        .rename(columns={"size": "frames"})
    )
    total = occ.groupby(["system", "replicate"])["frames"].transform("sum")
    occ["occupancy"] = occ["frames"] / total
    occ["centroid_M120_distance"] = occ["basin"].map(lambda b: centroids_map[b][0])
    occ["centroid_DRY_PIF_distance"] = occ["basin"].map(lambda b: centroids_map[b][1])
    return fes, occ.round(5)


def parse_path_output(path: Path, system: str) -> pd.DataFrame:
    rows = []
    if not path.exists():
        return pd.DataFrame()
    pattern = re.compile(r"Path:\s*\[(.*?)\],\s*Total Weight:\s*([0-9.eE+-]+)")
    with path.open(errors="ignore") as handle:
        for rank, line in enumerate(handle, start=1):
            m = pattern.search(line)
            if not m:
                continue
            raw_path = re.sub(r"np\.int64\((\d+)\)", r"\1", m.group(1))
            residues = [int(x) for x in re.findall(r"\d+", raw_path)]
            rows.append(
                {
                    "system": system,
                    "rank": rank,
                    "path_residues": "-".join(map(str, residues)),
                    "path_length": len(residues),
                    "total_weight": float(m.group(2)),
                    "contains_Y258_equiv": any(r in residues for r in [225, 227]),
                    "contains_M120_equiv": any(r in residues for r in [87, 89]),
                    "contains_Y222_equiv": any(r in residues for r in [189, 191]),
                }
            )
    return pd.DataFrame(rows)


def collect_mdpath_summary() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    md_roots = {
        "BM213": ROOT / "gpcr/bm213/mdpath/mdpath_out",
        "C5apep": ROOT / "gpcr/c5apep/mdpath/mdpath_out",
    }
    path_rows = []
    nmi_rows = []
    ci_rows = []
    for system, folder in md_roots.items():
        path_rows.append(parse_path_output(folder / "output.txt", system))
        nmi = pd.read_csv(folder / "nmi_df.csv")
        nmi["system"] = system
        nmi_rows.append(nmi)
        ci_path = folder / "path_confidence_intervals.txt"
        if ci_path.exists():
            with ci_path.open(errors="ignore") as handle:
                for line in handle:
                    m = re.search(r"^(.*?):\s*Mean=([-0-9.eE+]+),\s*2.5%=([-0-9.eE+]+),\s*97.5%=([-0-9.eE+]+)", line.strip())
                    if m:
                        ci_rows.append(
                            {
                                "system": system,
                                "path": m.group(1),
                                "mean": float(m.group(2)),
                                "ci_2p5": float(m.group(3)),
                                "ci_97p5": float(m.group(4)),
                            }
                        )
    paths = pd.concat(path_rows, ignore_index=True) if path_rows else pd.DataFrame()
    nmi_all = pd.concat(nmi_rows, ignore_index=True) if nmi_rows else pd.DataFrame()
    ci = pd.DataFrame(ci_rows)

    nmi_summary = (
        nmi_all.groupby("system", as_index=False)["MI Difference"]
        .agg(n_pairs="count", mean_mi_diff="mean", median_mi_diff="median", max_mi_diff="max")
        .round(6)
    )
    path_summary = (
        paths.groupby("system", as_index=False)
        .agg(top_weight=("total_weight", "max"), mean_top20_weight=("total_weight", lambda x: x.head(20).mean()), n_paths=("rank", "count"), mean_path_length=("path_length", "mean"))
        .round(5)
    )
    summary = path_summary.merge(nmi_summary, on="system", how="outer")
    return summary, paths, ci


def copy_residue_mapping() -> None:
    src = ROOT / "article/structure/04_reference_mapping/residue_mapping.csv"
    if src.exists():
        shutil.copyfile(src, TABLE_DIR / "residue_mapping.csv")


def main() -> None:
    ensure_dirs()
    warnings: list[str] = []
    write_placeholder_readme()

    hbonds = pd.concat([collect_metric(metric, warnings) for metric in HBOND_METRICS], ignore_index=True)
    microswitch = pd.concat([collect_metric(metric, warnings) for metric in MICROSWITCH_METRICS], ignore_index=True)
    rmsd = pd.concat([collect_metric(metric, warnings) for metric in RMSD_METRICS], ignore_index=True)

    hbonds.to_csv(DATA_DIR / "hbond_timeseries.csv", index=False)
    microswitch.to_csv(DATA_DIR / "microswitch_timeseries.csv", index=False)
    rmsd.to_csv(DATA_DIR / "rmsd_timeseries.csv", index=False)

    summarize_series(hbonds, ["system", "replicate", "metric"]).to_csv(TABLE_DIR / "hbond_summary.csv", index=False)
    summarize_series(microswitch, ["system", "replicate", "metric"]).to_csv(TABLE_DIR / "microswitch_summary.csv", index=False)

    inventory = collect_inventory([hbonds, microswitch, rmsd])
    inventory.to_csv(TABLE_DIR / "trajectory_inventory.csv", index=False)

    qc = collect_qc()
    qc.to_csv(TABLE_DIR / "production_qc.csv", index=False)

    fes, basin_summary = make_fes_points(microswitch, warnings)
    fes.to_csv(DATA_DIR / "fes_points.csv", index=False)
    basin_summary.to_csv(TABLE_DIR / "fes_basin_summary.csv", index=False)

    md_summary, md_paths, md_ci = collect_mdpath_summary()
    md_summary.to_csv(TABLE_DIR / "mdpath_summary.csv", index=False)
    md_paths.to_csv(DATA_DIR / "mdpath_top_paths.csv", index=False)
    md_ci.to_csv(DATA_DIR / "mdpath_path_confidence.csv", index=False)

    copy_residue_mapping()
    if not (TABLE_DIR / "residue_mapping.csv").exists():
        warnings.append("Residue mapping table not copied because article/structure mapping is absent.")

    warnings.append("PNG-only export policy active: PDF/SVG/TIFF were intentionally not generated.")
    warnings.append("Water-network, apo, C5a, PMX53, c5a_G, and arrestin docking inputs were intentionally excluded from this figure package.")
    note_warning(warnings)


if __name__ == "__main__":
    main()
