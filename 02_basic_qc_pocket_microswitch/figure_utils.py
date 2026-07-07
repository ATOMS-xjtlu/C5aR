from __future__ import annotations

import csv
import json
import math
import os
import re
from datetime import datetime
from pathlib import Path
from typing import Iterable

ROOT = Path(__file__).resolve().parents[3]
FIG_ROOT = ROOT / "article" / "figure"
ANALYSIS_ROOT = ROOT / "article" / "analysis" / "basic"
DATA_DIR = FIG_ROOT / "data"
TABLE_DIR = FIG_ROOT / "tables"
MAIN_DIR = FIG_ROOT / "main"
SUPP_DIR = FIG_ROOT / "supplement"
LOG_DIR = FIG_ROOT / "logs"
PLACEHOLDER_DIR = FIG_ROOT / "render_placeholders"
BASIC_MAIN_DIR = ANALYSIS_ROOT / "main"
BASIC_SUPP_DIR = ANALYSIS_ROOT / "supplement"

os.environ.setdefault("MPLCONFIGDIR", str(LOG_DIR / "mplconfig"))

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image
from scipy.cluster.vq import kmeans2, vq
from scipy.ndimage import gaussian_filter

BM_COLOR = "#00A6A6"
PEP_COLOR = "#EE7733"
NEUTRAL = "#6B7280"
ACCENT = "#009988"
WARN = "#CC3311"
SYSTEM_COLORS = {"BM213": BM_COLOR, "C5apep": PEP_COLOR}

SYSTEMS = [
    {"system": "BM213", "replicate": "rep1", "dir": ROOT / "gpcr/bm213", "prefix": "bm1"},
    {"system": "BM213", "replicate": "rep2", "dir": ROOT / "gpcr/bm213/second", "prefix": "bm2"},
    {"system": "BM213", "replicate": "rep3", "dir": ROOT / "gpcr/bm213/third", "prefix": "bm3"},
    {"system": "C5apep", "replicate": "rep1", "dir": ROOT / "gpcr/c5apep", "prefix": "pep1"},
    {"system": "C5apep", "replicate": "rep2", "dir": ROOT / "gpcr/c5apep/second", "prefix": "pep2"},
    {"system": "C5apep", "replicate": "rep3", "dir": ROOT / "gpcr/c5apep/third", "prefix": "pep3"},
]

METRIC_FILES = {
    ("BM213", "rep1", "Y258_hbond"): "258hb_bm1.dat",
    ("BM213", "rep2", "Y258_hbond"): "258hb_bm2.dat",
    ("BM213", "rep3", "Y258_hbond"): "258hb_bm3.dat",
    ("C5apep", "rep1", "Y258_hbond"): "258hb_pep1.dat",
    ("C5apep", "rep2", "Y258_hbond"): "258hb_pep2.dat",
    ("C5apep", "rep3", "Y258_hbond"): "258hb_pep3.dat",
    ("BM213", "rep1", "S171_hbond"): "171hb_bm1.dat",
    ("BM213", "rep2", "S171_hbond"): "171hb_bm2.dat",
    ("BM213", "rep3", "S171_hbond"): "171hb_bm3.dat",
    ("C5apep", "rep1", "S171_hbond"): "171hb_pep1.dat",
    ("C5apep", "rep2", "S171_hbond"): "171hb_pep2.dat",
    ("C5apep", "rep3", "S171_hbond"): "171hb_pep3.dat",
    ("BM213", "rep1", "M120_distance"): "M120ori_dist_bm1.dat",
    ("BM213", "rep2", "M120_distance"): "M120ori_dist_bm2.dat",
    ("BM213", "rep3", "M120_distance"): "M120ori_dist_bm3.dat",
    ("C5apep", "rep1", "M120_distance"): "M120ori_dist_pep1",
    ("C5apep", "rep2", "M120_distance"): "M120ori_dist_pep2.dat",
    ("C5apep", "rep3", "M120_distance"): "M120ori_dist_pep3.dat",
    ("BM213", "rep1", "W255_chi2"): "W255ori_bm1_X2.dat",
    ("BM213", "rep2", "W255_chi2"): "W255ori_bm2_X2.dat",
    ("BM213", "rep3", "W255_chi2"): "W255ori_bm3_X2.dat",
    ("C5apep", "rep1", "W255_chi2"): "W255ori_c5apep1_X2.dat",
    ("C5apep", "rep2", "W255_chi2"): "W255ori_c5apep2_X2.dat",
    ("C5apep", "rep3", "W255_chi2"): "W255ori_c5apep3_X2.dat",
    ("BM213", "rep1", "R350_chi2"): "R3.50_bm213_1.dat",
    ("BM213", "rep2", "R350_chi2"): "R3.50_bm213_2.dat",
    ("BM213", "rep3", "R350_chi2"): "R3.50_bm213_3.dat",
    ("C5apep", "rep1", "R350_chi2"): "R3.50_pep1_X2.dat",
    ("C5apep", "rep2", "R350_chi2"): "R3.50_pep2_X2.dat",
    ("C5apep", "rep3", "R350_chi2"): "R3.50_pep3_X2.dat",
    ("BM213", "rep1", "TM62_IC_distance"): "TM62_IC_bm1.dat",
    ("BM213", "rep2", "TM62_IC_distance"): "TM62_IC_bm2.dat",
    ("BM213", "rep3", "TM62_IC_distance"): "TM62_IC_bm3.dat",
    ("C5apep", "rep1", "TM62_IC_distance"): "TM62_IC_pep1.dat",
    ("C5apep", "rep2", "TM62_IC_distance"): "TM62_IC_pep2.dat",
    ("C5apep", "rep3", "TM62_IC_distance"): "TM62_IC_pep3.dat",
    ("BM213", "rep1", "TM72_IC_distance"): "TM72_IC_bm1.dat",
    ("BM213", "rep2", "TM72_IC_distance"): "TM72_IC_bm2.dat",
    ("BM213", "rep3", "TM72_IC_distance"): "TM72_IC_bm3.dat",
    ("C5apep", "rep1", "TM72_IC_distance"): "TM72_IC_pep1.dat",
    ("C5apep", "rep2", "TM72_IC_distance"): "TM72_IC_pep2.dat",
    ("C5apep", "rep3", "TM72_IC_distance"): "TM72_IC_pep3.dat",
    ("BM213", "rep1", "DRY_PIF_distance"): "dist-189-218.dat",
    ("BM213", "rep2", "DRY_PIF_distance"): "DRY_PIF_bm2.dat",
    ("BM213", "rep3", "DRY_PIF_distance"): "DRY_PIF_bm3.dat",
    ("C5apep", "rep1", "DRY_PIF_distance"): "DRY_PIF_pep1.dat",
    ("C5apep", "rep2", "DRY_PIF_distance"): "DRY_PIF_pep2.dat",
    ("C5apep", "rep3", "DRY_PIF_distance"): "DRY_PIF_pep3.dat",
    ("BM213", "rep1", "RMSD"): "rmsd-bm1.dat",
    ("BM213", "rep2", "RMSD"): "rmsd-bm2.dat",
    ("BM213", "rep3", "RMSD"): "rmsd-bm3.dat",
    ("C5apep", "rep1", "RMSD"): "rmsd-pep1.dat",
    ("C5apep", "rep2", "RMSD"): "rmsd-pep2.dat",
    ("C5apep", "rep3", "RMSD"): "rmsd-pep3.dat",
}

METRIC_LABELS = {
    "Y258_hbond": "Y258-ligand H-bond",
    "S171_hbond": "S171-ligand H-bond",
    "M120_distance": "M120-D82 distance (Å)",
    "W255_chi2": "W255 χ2 (deg)",
    "R350_chi2": "R3.50 χ2 (deg)",
    "R3.50_chi2": "R3.50 χ2 (deg)",
    "TM62_IC_distance": "TM6-TM2 IC distance (Å)",
    "TM72_IC_distance": "TM7-TM2 IC distance (Å)",
    "DRY_PIF_distance": "Y222-F251 distance (Å)",
    "RMSD": "RMSD (Å)",
}


def ensure_dirs() -> None:
    for path in [
        DATA_DIR,
        TABLE_DIR,
        MAIN_DIR,
        SUPP_DIR,
        LOG_DIR,
        PLACEHOLDER_DIR,
        BASIC_MAIN_DIR,
        BASIC_SUPP_DIR,
        LOG_DIR / "mplconfig",
    ]:
        path.mkdir(parents=True, exist_ok=True)


def apply_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "Liberation Sans", "sans-serif"],
            "font.size": 12,
            "axes.labelsize": 16,
            "axes.titlesize": 15,
            "xtick.labelsize": 12,
            "ytick.labelsize": 12,
            "legend.fontsize": 11,
            "axes.spines.right": True,
            "axes.spines.top": True,
            "axes.linewidth": 1.0,
            "axes.grid": False,
            "legend.frameon": False,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.dpi": 600,
        }
    )


def save_png(fig: plt.Figure, path: Path, dpi: int = 600) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def save_pngs(fig: plt.Figure, paths: Iterable[Path], dpi: int = 600) -> None:
    for path in paths:
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.08, 1.05, label, transform=ax.transAxes, ha="left", va="bottom", fontweight="bold", fontsize=12)


def style_axis(ax: plt.Axes) -> None:
    ax.tick_params(which="both", direction="out", top=False, right=False, bottom=True, left=True, length=4.5, width=1.0, pad=3)
    for spine in ax.spines.values():
        spine.set_linewidth(1.0)


def metric_path(system: str, replicate: str, metric: str) -> Path:
    info = next(x for x in SYSTEMS if x["system"] == system and x["replicate"] == replicate)
    name = METRIC_FILES[(system, replicate, metric)]
    return info["dir"] / name


def read_dat(path: Path) -> pd.DataFrame:
    rows = []
    with path.open(errors="ignore") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            parts = stripped.split()
            if len(parts) < 2:
                continue
            try:
                rows.append((int(float(parts[0])), float(parts[1])))
            except ValueError:
                continue
    return pd.DataFrame(rows, columns=["frame", "value"])


def add_time_ns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if out.empty:
        out["time_ns"] = []
        return out
    max_frame = max(int(out["frame"].max()), 1)
    if max_frame == 1:
        out["time_ns"] = 0.0
    else:
        out["time_ns"] = (out["frame"] - 1) * (2000.0 / (max_frame - 1))
    return out


def downsample(df: pd.DataFrame, max_points: int = 2500) -> pd.DataFrame:
    if len(df) <= max_points:
        return df
    idx = np.linspace(0, len(df) - 1, max_points).astype(int)
    return df.iloc[idx]


def rolling_mean(values: Iterable[float], window: int = 250) -> np.ndarray:
    arr = np.asarray(list(values), dtype=float)
    if len(arr) < window:
        return arr
    kernel = np.ones(window) / window
    return np.convolve(arr, kernel, mode="same")


def summarize_series(df: pd.DataFrame, group_cols: list[str], value_col: str = "value") -> pd.DataFrame:
    return (
        df.groupby(group_cols, as_index=False)[value_col]
        .agg(frames="count", mean="mean", sd="std", median="median", q25=lambda x: x.quantile(0.25), q75=lambda x: x.quantile(0.75), minimum="min", maximum="max")
        .round(5)
    )


def parse_production_log(path: Path) -> dict:
    record = {
        "log_path": str(path.relative_to(ROOT)),
        "exists": path.exists(),
        "nstlim": np.nan,
        "dt_ps": np.nan,
        "temp0_k": np.nan,
        "ntwx": np.nan,
        "last_nstep": np.nan,
        "last_time_ps": np.nan,
        "mean_temp_k": np.nan,
        "mean_density": np.nan,
        "error_or_warning_count": 0,
    }
    if not path.exists():
        return record
    temps = []
    densities = []
    first_lines = []
    with path.open(errors="ignore") as handle:
        for _ in range(1400):
            line = handle.readline()
            if not line:
                break
            first_lines.append(line)
    text_head = "".join(first_lines)
    for key, regex, cast in [
        ("nstlim", r"nstlim\s*=\s*([0-9]+)", int),
        ("dt_ps", r"dt\s*=\s*([0-9.]+)", float),
        ("temp0_k", r"temp0\s*=\s*([0-9.]+)", float),
        ("ntwx", r"ntwx\s*=\s*([0-9]+)", int),
    ]:
        m = re.search(regex, text_head)
        if m:
            record[key] = cast(m.group(1))

    with path.open("rb") as handle:
        handle.seek(0, os.SEEK_END)
        size = handle.tell()
        handle.seek(max(size - 2_000_000, 0))
        tail_text = handle.read().decode(errors="ignore")
    text_sample = text_head + "\n" + tail_text
    for line in text_sample.splitlines():
        if "NSTEP =" in line and "TEMP(K)" in line:
            m = re.search(r"NSTEP\s*=\s*([0-9]+).*TIME\(PS\)\s*=\s*([-0-9.]+).*TEMP\(K\)\s*=\s*([-0-9.]+)", line)
            if m:
                record["last_nstep"] = int(m.group(1))
                record["last_time_ps"] = float(m.group(2))
                temps.append(float(m.group(3)))
        if "Density" in line:
            m = re.search(r"Density\s*=\s*([-0-9.]+)", line)
            if m:
                densities.append(float(m.group(1)))
        if re.search(r"\b(ERROR|WARNING|NaN|nan)\b", line):
            record["error_or_warning_count"] += 1
    if temps:
        record["mean_temp_k"] = float(np.mean(temps))
    if densities:
        record["mean_density"] = float(np.mean(densities))
    return record


def load_required_csv(name: str) -> pd.DataFrame:
    path = DATA_DIR / name
    if not path.exists():
        raise FileNotFoundError(f"Run collect_data.py first; missing {path}")
    return pd.read_csv(path)


def draw_placeholder(ax: plt.Axes, title: str, lines: list[str]) -> None:
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("#111827")
        spine.set_linewidth(1.0)
    ax.set_facecolor("#F7F7F7")
    ax.text(0.05, 0.88, title, transform=ax.transAxes, fontweight="bold", fontsize=10, va="top")
    ax.text(0.05, 0.72, "\n".join(lines), transform=ax.transAxes, fontsize=8, va="top", linespacing=1.45)


def grouped_bar_points(ax: plt.Axes, summary: pd.DataFrame, metric: str, ylabel: str, ylim: tuple[float, float] | None = None) -> None:
    sub = summary[summary["metric"] == metric].copy()
    systems = ["BM213", "C5apep"]
    x = np.arange(len(systems))
    means = [sub[sub["system"] == s]["mean"].mean() for s in systems]
    sds = [sub[sub["system"] == s]["mean"].std(ddof=1) for s in systems]
    ax.bar(x, means, yerr=sds, color=[SYSTEM_COLORS[s] for s in systems], edgecolor="black", width=0.62, linewidth=0.8, capsize=3)
    for i, system in enumerate(systems):
        vals = sub[sub["system"] == system].sort_values("replicate")["mean"].to_numpy()
        jitter = np.linspace(-0.12, 0.12, len(vals)) if len(vals) > 1 else [0]
        ax.scatter(i + jitter, vals, s=28, color="white", edgecolor="black", zorder=4)
    ax.set_xticks(x)
    ax.set_xticklabels(systems)
    ax.set_ylabel(ylabel)
    if ylim:
        ax.set_ylim(*ylim)
    style_axis(ax)


def plot_kde_or_hist(ax: plt.Axes, values: np.ndarray, color: str, label: str, bins: int = 80, density: bool = True) -> None:
    clean = np.asarray(values, dtype=float)
    clean = clean[np.isfinite(clean)]
    if len(clean) == 0:
        return
    counts, edges = np.histogram(clean, bins=bins, density=density)
    centers = (edges[:-1] + edges[1:]) / 2
    if len(counts) > 5:
        counts = gaussian_filter(counts.astype(float), sigma=1.2)
    ax.plot(centers, counts, color=color, lw=1.8, label=label)


def make_free_energy(x: np.ndarray, y: np.ndarray, bins: int = 55, ranges: tuple[tuple[float, float], tuple[float, float]] | None = None):
    hist, xedges, yedges = np.histogram2d(x, y, bins=bins, range=ranges)
    hist = gaussian_filter(hist, sigma=1.0)
    prob = hist / np.nansum(hist)
    with np.errstate(divide="ignore", invalid="ignore"):
        fes = -0.593 * np.log(prob / np.nanmax(prob))
    fes[~np.isfinite(fes)] = np.nan
    return fes.T, xedges, yedges


def image_panel(ax: plt.Axes, path: Path, title: str) -> None:
    ax.set_xticks([])
    ax.set_yticks([])
    if path.exists():
        img = Image.open(path)
        ax.imshow(img)
    else:
        draw_placeholder(ax, title, [f"Missing image: {path.name}"])
    ax.set_title(title, fontsize=9, pad=3)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.8)


def write_placeholder_readme() -> None:
    text = """# Render placeholders

This figure set is PNG-only and keeps structural render panels as placeholders.

Required manual render replacements:

1. Pocket close-up: use `article/structure/01_pocket_closeup/*.pml`.
2. Representative basin/state structures: use `article/structure/02_representative_states/load_representative_states.pml`.
3. BM213/C5apep MDPath 3D pathway view: use `article/structure/03_mdpath_3d/load_*_mdpath_*.pml`.

The quantitative PNG panels in `article/figure/main` and `article/figure/supplement` can be regenerated with:

`python3 article/figure/scripts/run_all.py`

No PDF, SVG, or TIFF files are generated by this package.
"""
    (PLACEHOLDER_DIR / "README.md").write_text(text)


def expected_outputs() -> list[tuple[str, Path, str, str]]:
    return [
        ("Fig1_overview", MAIN_DIR / "Fig1_overview.png", "plot_fig1_overview.py", "tables plus schematic placeholders"),
        ("Fig2_pocket_contacts", MAIN_DIR / "Fig2_pocket_contacts.png", "plot_fig2_pocket_contacts.py", "hbond_timeseries.csv, hbond_summary.csv"),
        ("Fig2_pocket_contacts_basic", BASIC_MAIN_DIR / "Fig2_pocket_contacts.png", "plot_fig2_pocket_contacts.py", "hbond tables plus IWI BSA"),
        ("Fig2_iwi_bsa", BASIC_MAIN_DIR / "Fig2_iwi_bsa.png", "plot_fig2_pocket_contacts.py", "IWI pocket SASA component files"),
        ("Fig3_microswitch", MAIN_DIR / "Fig3_microswitch.png", "plot_fig3_microswitch.py", "microswitch_timeseries.csv, microswitch_summary.csv"),
        ("Fig4_fes_projection", MAIN_DIR / "Fig4_fes_projection.png", "plot_fig4_fes_projection.py", "fes_points.csv, fes_basin_summary.csv"),
        ("Fig5_mdpath_network", MAIN_DIR / "Fig5_mdpath_network.png", "plot_fig5_mdpath_network.py", "mdpath_summary.csv plus diagnostic PNGs"),
        ("FigS1_qc", SUPP_DIR / "FigS1_qc.png", "plot_figS1_qc.py", "rmsd_timeseries.csv, production_qc.csv"),
        ("C5AR_rmsf_regions", BASIC_SUPP_DIR / "C5AR_rmsf_regions.png", "plot_c5ar_rmsf_regions.py", "available RMSF windows from apo/C5a/C5a_G/BM213"),
        ("FigS2_hbond_timeseries", SUPP_DIR / "FigS2_hbond_timeseries.png", "plot_figS2_hbond_timeseries.py", "hbond_timeseries.csv"),
        ("FigS3_microswitch_replicates", SUPP_DIR / "FigS3_microswitch_replicates.png", "plot_figS3_microswitch_replicates.py", "microswitch_timeseries.csv"),
        ("FigS4_fes_basin_definition", SUPP_DIR / "FigS4_fes_basin_definition.png", "plot_figS4_fes_basin_definition.py", "fes_points.csv, fes_basin_summary.csv"),
        ("FigS5_mdpath_diagnostics", SUPP_DIR / "FigS5_mdpath_diagnostics.png", "plot_figS5_mdpath_diagnostics.py", "MDPath graph/clustered_paths PNGs, nmi_df.csv, path CI"),
    ]


def write_figure_manifest() -> None:
    rows = []
    now = datetime.now().isoformat(timespec="seconds")
    for figure, out, script, inputs in expected_outputs():
        rows.append(
            {
                "figure": figure,
                "output_path": str(out.relative_to(ROOT)),
                "script": f"article/figure/scripts/{script}",
                "input_summary": inputs,
                "format": "png",
                "generated_at": now,
                "exists": out.exists(),
                "size_bytes": out.stat().st_size if out.exists() else 0,
            }
        )
    pd.DataFrame(rows).to_csv(TABLE_DIR / "figure_manifest.csv", index=False)


def note_warning(lines: list[str]) -> None:
    if not lines:
        lines = ["No data warnings generated."]
    (LOG_DIR / "data_warnings.md").write_text("# Data warnings\n\n" + "\n".join(f"- {line}" for line in lines) + "\n")
