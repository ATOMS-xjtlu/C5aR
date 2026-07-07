from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path("article/figure/logs/mplconfig").resolve()))

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd
import MDAnalysis as mda
from Bio.PDB.Atom import Atom
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Residue import Residue
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.Structure import Structure

from figure_utils import ANALYSIS_ROOT, DATA_DIR, SYSTEM_COLORS, TABLE_DIR, apply_style, grouped_bar_points, plot_kde_or_hist, save_png, style_axis


OUT_DIRS = {
    "data": ANALYSIS_ROOT / "data",
    "qc": ANALYSIS_ROOT / "qc",
    "pocket": ANALYSIS_ROOT / "pocket",
    "microswitch": ANALYSIS_ROOT / "microswitch",
    "flexibility": ANALYSIS_ROOT / "flexibility",
}

RMSD_WINDOW_NS = 2000.0

APO_COLOR = "#6B7280"
SYSTEM_ORDER = ["apo", "BM213", "C5apep"]
DISPLAY = {"apo": "apo", "BM213": "BM213", "C5apep": "C5apep"}
DISPLAY_COLORS = {"apo": APO_COLOR, "BM213": SYSTEM_COLORS["BM213"], "C5apep": SYSTEM_COLORS["C5apep"]}


@dataclass(frozen=True)
class TrajectorySpec:
    system: str
    replicate: str
    topology: Path
    trajectory: Path
    c5ar_offset: int
    total_ns: float


TRAJECTORIES = [
    TrajectorySpec("apo", "rep1", Path("gpcr/apo/md-initial.pdb"), Path("article/analysis/alloviz/inputs/downsampled/formal/apo_rep1_stride10.xtc"), 21, 2000.0),
    TrajectorySpec("apo", "rep2", Path("gpcr/apo/second/md-initial.pdb"), Path("article/analysis/alloviz/inputs/downsampled/formal/apo_rep2_stride10.xtc"), 21, 2000.0),
    TrajectorySpec("apo", "rep3", Path("gpcr/apo/third/md-initial.pdb"), Path("article/analysis/alloviz/inputs/downsampled/formal/apo_rep3_stride10.xtc"), 21, 2000.0),
    TrajectorySpec("BM213", "rep1", Path("gpcr/bm213/md-initial.pdb"), Path("article/analysis/alloviz/inputs/downsampled/formal/BM213_rep1_stride10.xtc"), 33, 2000.0),
    TrajectorySpec("BM213", "rep2", Path("gpcr/bm213/second/md-initial.pdb"), Path("article/analysis/alloviz/inputs/downsampled/formal/BM213_rep2_stride10.xtc"), 33, 2000.0),
    TrajectorySpec("BM213", "rep3", Path("gpcr/bm213/third/md-initial.pdb"), Path("article/analysis/alloviz/inputs/downsampled/formal/BM213_rep3_stride10.xtc"), 33, 2000.0),
    TrajectorySpec("C5apep", "rep1", Path("gpcr/c5apep/md-initial.pdb"), Path("article/analysis/alloviz/inputs/downsampled/formal/C5apep_rep1_stride10.xtc"), 31, 2000.0),
    TrajectorySpec("C5apep", "rep2", Path("gpcr/c5apep/second/md-initial.pdb"), Path("article/analysis/alloviz/inputs/downsampled/formal/C5apep_rep2_stride10.xtc"), 31, 2000.0),
    TrajectorySpec("C5apep", "rep3", Path("gpcr/c5apep/third/md-initial.pdb"), Path("article/analysis/alloviz/inputs/downsampled/formal/C5apep_rep3_stride10.xtc"), 31, 3000.0),
]

IWI_SPECS = {
    "BM213": {
        "pocket_resids": [58, 69, 83],
        "ligand_resids": list(range(282, 290)),
    },
    "C5apep": {
        "pocket_resids": [60, 71, 85],
        "ligand_resids": [284],
    },
}

MICRO_METRICS = [
    ("M120_distance", "m120_d82_by_system.png", "M120-D82 distance (Å)", "M120-D82"),
    ("DRY_PIF_distance", "y222_f251_by_system.png", "Y222-F251 distance (Å)", "Y222-F251"),
    ("W255_chi2", "w255_chi2_by_system.png", "W255 χ2 (deg)", "W255 χ2"),
    ("R350_chi2", "r350_chi2_by_system.png", "R134 χ2 (deg)", "R134 χ2"),
    ("TM62_IC_distance", "tm6_tm2_ic_by_system.png", "TM6-TM2 IC distance (Å)", "TM6-TM2"),
    ("TM72_IC_distance", "tm7_tm2_ic_by_system.png", "TM7-TM2 IC distance (Å)", "TM7-TM2"),
]

TM_REGION_BANDS = [
    ("TM1", 21, 49, "#D9F0F0", "tm"),
    ("ECL1", 50, 57, "#E8ECEF", "loop"),
    ("TM2", 58, 84, "#D9F0F0", "tm"),
    ("ICL1", 85, 90, "#E8ECEF", "loop"),
    ("TM3", 91, 118, "#D9F0F0", "tm"),
    ("ECL2", 119, 139, "#E8ECEF", "loop"),
    ("TM4", 140, 170, "#D9F0F0", "tm"),
    ("ICL2", 171, 188, "#E8ECEF", "loop"),
    ("TM5", 189, 225, "#D9F0F0", "tm"),
    ("ECL3", 226, 239, "#E8ECEF", "loop"),
    ("TM6", 240, 266, "#D9F0F0", "tm"),
    ("ICL3", 267, 276, "#E8ECEF", "loop"),
    ("TM7", 277, 300, "#D9F0F0", "tm"),
    ("H8", 301, 314, "#CFEDED", "tm"),
]

POCKET_SINGLE_FIGSIZE = (4.5, 3.9)
POCKET_WIDE_FIGSIZE = (6.2, 4.0)
POCKET_LABEL_SIZE = 13
POCKET_TICK_SIZE = 10.5
POCKET_HBOND_YLIM = (0, 1.8)
POCKET_HBOND_YTICKS = np.arange(0.0, 1.51, 0.5)


def ensure_output_dirs() -> None:
    for path in OUT_DIRS.values():
        path.mkdir(parents=True, exist_ok=True)


def kabsch_align(mobile: np.ndarray, reference: np.ndarray) -> np.ndarray:
    mob_center = mobile.mean(axis=0)
    ref_center = reference.mean(axis=0)
    p = mobile - mob_center
    q = reference - ref_center
    cov = p.T @ q
    v, _s, wt = np.linalg.svd(cov)
    if np.linalg.det(v @ wt) < 0:
        v[:, -1] *= -1
    rot = v @ wt
    return p @ rot + ref_center


def calc_rmsd_rmsf(spec: TrajectorySpec) -> tuple[pd.DataFrame, pd.DataFrame]:
    u = mda.Universe(str(spec.topology), str(spec.trajectory))
    ca = u.select_atoms("protein and name CA")
    if len(ca) == 0:
        raise ValueError(f"No protein CA atoms for {spec.system} {spec.replicate}")
    u.trajectory[0]
    reference = ca.positions.astype(np.float64).copy()
    mean = np.zeros_like(reference)
    m2 = np.zeros_like(reference)
    rmsd_rows = []
    n_frames = len(u.trajectory)
    for idx, ts in enumerate(u.trajectory):
        aligned = kabsch_align(ca.positions.astype(np.float64), reference)
        diff_ref = aligned - reference
        rmsd = float(np.sqrt(np.mean(np.sum(diff_ref * diff_ref, axis=1))))
        time_ns = (idx / max(n_frames - 1, 1)) * spec.total_ns
        rmsd_rows.append((spec.system, spec.replicate, idx + 1, time_ns, rmsd))
        count = idx + 1
        delta = aligned - mean
        mean += delta / count
        m2 += delta * (aligned - mean)
    variance = np.sum(m2, axis=1) / max(n_frames - 1, 1)
    rmsf = np.sqrt(variance)
    rmsf_df = pd.DataFrame(
        {
            "system": spec.system,
            "replicate": spec.replicate,
            "local_resid": ca.resids,
            "c5ar_resid": ca.resids + spec.c5ar_offset,
            "resname": ca.resnames,
            "rmsf_a": rmsf,
        }
    )
    rmsd_df = pd.DataFrame(rmsd_rows, columns=["system", "replicate", "frame", "time_ns", "rmsd_a"])
    return rmsd_df, rmsf_df


def build_ca_dynamics() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    rmsd_path = OUT_DIRS["data"] / "ca_rmsd_timeseries.csv"
    rmsf_path = OUT_DIRS["data"] / "ca_rmsf_by_replicate.csv"
    rmsf_summary_path = OUT_DIRS["data"] / "ca_rmsf_system_summary.csv"
    if rmsd_path.exists() and rmsf_path.exists() and rmsf_summary_path.exists():
        return pd.read_csv(rmsd_path), pd.read_csv(rmsf_path), pd.read_csv(rmsf_summary_path)
    rmsd_blocks = []
    rmsf_blocks = []
    for spec in TRAJECTORIES:
        rmsd, rmsf = calc_rmsd_rmsf(spec)
        rmsd_blocks.append(rmsd)
        rmsf_blocks.append(rmsf)
    rmsd_df = pd.concat(rmsd_blocks, ignore_index=True)
    rmsf_df = pd.concat(rmsf_blocks, ignore_index=True)
    rmsf_summary = (
        rmsf_df.groupby(["system", "c5ar_resid"], as_index=False)
        .agg(mean_rmsf_a=("rmsf_a", "mean"), sd_rmsf_a=("rmsf_a", "std"), n_reps=("replicate", "nunique"))
        .sort_values(["system", "c5ar_resid"])
    )
    rmsd_df.to_csv(rmsd_path, index=False)
    rmsf_df.to_csv(rmsf_path, index=False)
    rmsf_summary.to_csv(rmsf_summary_path, index=False)
    return rmsd_df, rmsf_df, rmsf_summary


def smooth(values: np.ndarray, window: int = 7) -> np.ndarray:
    if len(values) < window:
        return values
    kernel = np.ones(window, dtype=float) / window
    return np.convolve(values, kernel, mode="same")


def plot_rmsd_timeseries(rmsd_df: pd.DataFrame) -> None:
    rmsd_df = rmsd_df[rmsd_df["time_ns"] <= RMSD_WINDOW_NS].copy()
    fig, ax = plt.subplots(figsize=(7.3, 5.0))
    for system in SYSTEM_ORDER:
        sub_system = rmsd_df[rmsd_df["system"] == system]
        for rep, sub in sub_system.groupby("replicate"):
            sub = sub.sort_values("time_ns")
            alpha = 0.55 if system != "apo" else 0.48
            ax.plot(sub["time_ns"], smooth(sub["rmsd_a"].to_numpy(), 9), color=DISPLAY_COLORS[system], lw=1.1, alpha=alpha, label="_nolegend_")
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("RMSD (Å)")
    ax.set_xlim(0, RMSD_WINDOW_NS)
    ax.set_ylim(0, 8)
    ax.set_yticks([0, 2, 4, 6, 8])
    legend_labels = {"apo": "apo", "BM213": "BM213", "C5apep": "C5Apep"}
    handles = [Line2D([0], [0], color=DISPLAY_COLORS[system], lw=2.0, label=legend_labels[system]) for system in SYSTEM_ORDER]
    ax.legend(handles=handles, fontsize=10.5, loc="upper right", frameon=True, framealpha=0.9, facecolor="white", edgecolor="none")
    style_axis(ax)
    save_png(fig, OUT_DIRS["qc"] / "FigS1_qc.png")


def plot_rmsd_distribution(rmsd_df: pd.DataFrame) -> None:
    rmsd_df = rmsd_df[rmsd_df["time_ns"] <= RMSD_WINDOW_NS].copy()
    fig, ax = plt.subplots(figsize=(5.6, 4.4))
    for system in SYSTEM_ORDER:
        vals = rmsd_df[rmsd_df["system"] == system]["rmsd_a"].to_numpy()
        plot_kde_or_hist(ax, vals, DISPLAY_COLORS[system], DISPLAY[system], bins=70)
    ax.set_xlabel("RMSD (Å)")
    ax.set_ylabel("Density")
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.legend(fontsize=8, loc="upper right")
    style_axis(ax)
    save_png(fig, OUT_DIRS["qc"] / "rmsd_distribution.png")


def plot_hbond_basic_panels() -> None:
    summary = pd.read_csv("article/figure/tables/hbond_summary.csv")
    specs = [
        ("Y258_hbond", "y258_ligand_hbond.png", "Mean H-bond count", "Y258-ligand", POCKET_HBOND_YLIM, POCKET_HBOND_YTICKS),
        ("S171_hbond", "s171_ligand_hbond.png", "Mean H-bond count", "S171-ligand", POCKET_HBOND_YLIM, POCKET_HBOND_YTICKS),
    ]
    for metric, filename, ylabel, xlabel, ylim, yticks in specs:
        fig, ax = plt.subplots(figsize=POCKET_SINGLE_FIGSIZE)
        grouped_bar_points(ax, summary, metric, ylabel, ylim)
        ax.set_ylabel(f"{xlabel} {ylabel}")
        style_pocket_axis(ax, y_bins=4)
        if yticks is not None:
            ax.set_yticks(yticks)
        save_png(fig, OUT_DIRS["pocket"] / filename)


def style_pocket_axis(ax: plt.Axes, y_bins: int = 5) -> None:
    ax.yaxis.set_major_locator(MaxNLocator(nbins=y_bins))
    ax.yaxis.label.set_size(POCKET_LABEL_SIZE)
    ax.xaxis.label.set_size(POCKET_LABEL_SIZE)
    ax.tick_params(axis="both", labelsize=POCKET_TICK_SIZE)


def infer_element(atom) -> str:
    element = getattr(atom, "element", "") or ""
    element = element.strip()
    if element:
        return element.upper()
    name = atom.name.strip()
    if len(name) >= 2 and name[:2].upper() in {"CL", "BR", "NA", "MG", "ZN", "FE", "CA"}:
        return name[:2].upper()
    return name[0].upper()


def bio_structure_for_atoms(name: str, coords: np.ndarray, atom_names: list[str], resids: list[int], resnames: list[str], elements: list[str]) -> Structure:
    structure = Structure(name)
    model = Model(0)
    chain = Chain("A")
    structure.add(model)
    model.add(chain)
    residues: dict[int, Residue] = {}
    for i, (coord, atom_name, resid, resname, element) in enumerate(zip(coords, atom_names, resids, resnames, elements), start=1):
        if resid not in residues:
            residue = Residue((" ", int(resid), " "), resname, " ")
            chain.add(residue)
            residues[resid] = residue
        fullname = atom_name[:4].rjust(4)
        residues[resid].add(Atom(atom_name[:4], coord.astype(float), 0.0, 1.0, " ", fullname, i, element=element))
    return structure


def sasa_for_selection(coords: np.ndarray, atom_names: list[str], resids: list[int], resnames: list[str], elements: list[str], sr: ShrakeRupley) -> float:
    structure = bio_structure_for_atoms("sasa", coords, atom_names, resids, resnames, elements)
    sr.compute(structure, level="S")
    return float(structure.sasa)


def calc_iwi_bsa_for_spec(spec: TrajectorySpec, max_samples: int = 450) -> pd.DataFrame:
    iwi = IWI_SPECS[spec.system]
    u = mda.Universe(str(spec.topology), str(spec.trajectory))
    pocket = u.select_atoms("resid " + " ".join(str(x) for x in iwi["pocket_resids"]))
    ligand = u.select_atoms("resid " + " ".join(str(x) for x in iwi["ligand_resids"]))
    if len(pocket) == 0 or len(ligand) == 0:
        raise ValueError(f"Missing IWI pocket or ligand selection for {spec.system} {spec.replicate}")
    complex_atoms = pocket + ligand
    frame_indices = np.linspace(0, len(u.trajectory) - 1, min(max_samples, len(u.trajectory))).astype(int)
    sr = ShrakeRupley(n_points=100)
    rows = []
    for idx in frame_indices:
        u.trajectory[int(idx)]
        time_ns = (idx / max(len(u.trajectory) - 1, 1)) * spec.total_ns
        sasa_pocket = sasa_for_selection(pocket.positions, list(pocket.names), list(pocket.resids), list(pocket.resnames), [infer_element(a) for a in pocket], sr)
        sasa_lig = sasa_for_selection(ligand.positions, list(ligand.names), list(ligand.resids), list(ligand.resnames), [infer_element(a) for a in ligand], sr)
        sasa_complex = sasa_for_selection(complex_atoms.positions, list(complex_atoms.names), list(complex_atoms.resids), list(complex_atoms.resnames), [infer_element(a) for a in complex_atoms], sr)
        buried = max(sasa_pocket + sasa_lig - sasa_complex, 0.0)
        rows.append((spec.system, spec.replicate, int(idx) + 1, time_ns, sasa_pocket, sasa_lig, sasa_complex, buried))
    return pd.DataFrame(rows, columns=["system", "replicate", "frame", "time_ns", "sasa_iwi", "sasa_ligand", "sasa_complex", "buried_sasa_a2"])


def build_iwi_bsa_replicates() -> pd.DataFrame:
    out_path = OUT_DIRS["data"] / "iwi_bsa_replicates.csv"
    summary_path = OUT_DIRS["data"] / "iwi_bsa_replicate_summary.csv"
    if out_path.exists() and summary_path.exists():
        return pd.read_csv(out_path)
    blocks = []
    for spec in TRAJECTORIES:
        if spec.system in IWI_SPECS:
            blocks.append(calc_iwi_bsa_for_spec(spec))
    df = pd.concat(blocks, ignore_index=True)
    summary = (
        df.groupby(["system", "replicate"], as_index=False)["buried_sasa_a2"]
        .agg(frames="count", mean="mean", sd="std", median="median", q25=lambda x: x.quantile(0.25), q75=lambda x: x.quantile(0.75), minimum="min", maximum="max")
        .round(5)
    )
    df.to_csv(out_path, index=False)
    summary.to_csv(summary_path, index=False)
    return df


def plot_iwi_bsa_replicates(bsa: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=POCKET_WIDE_FIGSIZE)
    positions = []
    vals = []
    colors = []
    labels = []
    pos = 1.0
    for system in ["BM213", "C5apep"]:
        for rep in ["rep1", "rep2", "rep3"]:
            sub = bsa[(bsa["system"] == system) & (bsa["replicate"] == rep)]["buried_sasa_a2"].dropna()
            vals.append(sub.to_numpy())
            positions.append(pos)
            colors.append(DISPLAY_COLORS[system])
            labels.append(f"{DISPLAY[system]}\n{rep}")
            pos += 1.0
        pos += 0.55
    bp = ax.boxplot(vals, positions=positions, patch_artist=True, widths=0.58, showfliers=False)
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.82)
        patch.set_edgecolor("black")
    for median in bp["medians"]:
        median.set_color("black")
        median.set_linewidth(1.1)
    ax.set_xticks(positions)
    ax.set_xticklabels(labels)
    ax.set_ylabel("IWI buried SASA (Å²)")
    style_axis(ax)
    style_pocket_axis(ax, y_bins=5)
    save_png(fig, OUT_DIRS["pocket"] / "iwi_bsa_replicates.png")


def plot_iwi_bsa_system_mean(bsa: pd.DataFrame) -> None:
    rep_means = (
        bsa.groupby(["system", "replicate"], as_index=False)["buried_sasa_a2"]
        .mean()
        .rename(columns={"buried_sasa_a2": "replicate_mean"})
    )
    fig, ax = plt.subplots(figsize=POCKET_SINGLE_FIGSIZE)
    systems = ["BM213", "C5apep"]
    positions = np.arange(1, len(systems) + 1)
    vals_by_system = [rep_means[rep_means["system"] == system]["replicate_mean"].to_numpy() for system in systems]
    bp = ax.boxplot(vals_by_system, positions=positions, patch_artist=True, widths=0.56, showfliers=False)
    for patch, system in zip(bp["boxes"], systems):
        patch.set_facecolor(DISPLAY_COLORS[system])
        patch.set_alpha(0.85)
        patch.set_edgecolor("black")
    for median in bp["medians"]:
        median.set_color("black")
        median.set_linewidth(1.1)
    ax.set_xticks(positions)
    ax.set_xticklabels([DISPLAY[system] for system in systems])
    ax.set_ylabel("Mean IWI buried SASA (Å²)")
    ax.set_ylim(0, 200)
    style_axis(ax)
    style_pocket_axis(ax, y_bins=4)
    ax.set_yticks([0, 50, 100, 150, 200])
    save_png(fig, OUT_DIRS["pocket"] / "iwi_bsa_mean.png")


def box_by_system(ax: plt.Axes, data: pd.DataFrame, metric: str, ylabel: str) -> None:
    vals = []
    for system in ["BM213", "C5apep"]:
        sub = data[(data["system"] == system) & (data["metric"] == metric)]["value"].dropna()
        vals.append(sub.sample(n=min(5000, len(sub)), random_state=4).to_numpy())
    bp = ax.boxplot(vals, patch_artist=True, widths=0.56, showfliers=False)
    for patch, system in zip(bp["boxes"], ["BM213", "C5apep"]):
        patch.set_facecolor(DISPLAY_COLORS[system])
        patch.set_alpha(0.85)
        patch.set_edgecolor("black")
    for median in bp["medians"]:
        median.set_color("black")
        median.set_linewidth(1.1)
    ax.set_xticks([1, 2])
    ax.set_xticklabels(["BM213", "C5apep"])
    ax.set_ylabel(ylabel)
    style_axis(ax)


def plot_microswitch_by_system() -> None:
    data = pd.read_csv("article/figure/data/microswitch_timeseries.csv")
    for metric, filename, ylabel, _title in MICRO_METRICS:
        fig, ax = plt.subplots(figsize=(4.3, 3.9))
        if "chi2" in metric:
            for system in ["BM213", "C5apep"]:
                sub = data[(data["system"] == system) & (data["metric"] == metric)]
                plot_kde_or_hist(ax, sub["value"].to_numpy(), DISPLAY_COLORS[system], system, bins=90)
            ax.set_xlabel(ylabel)
            ax.set_ylabel("Density")
            ax.legend(fontsize=9)
            style_axis(ax)
        else:
            box_by_system(ax, data, metric, ylabel)
        save_png(fig, OUT_DIRS["microswitch"] / filename)


def plot_microswitch_replicates() -> None:
    data = pd.read_csv("article/figure/data/microswitch_timeseries.csv")
    for metric, filename, ylabel, _title in MICRO_METRICS:
        fig, ax = plt.subplots(figsize=(6.0, 4.0))
        positions = []
        vals = []
        colors = []
        labels = []
        pos = 1.0
        for system in ["BM213", "C5apep"]:
            for rep in ["rep1", "rep2", "rep3"]:
                sub = data[(data["system"] == system) & (data["replicate"] == rep) & (data["metric"] == metric)]["value"].dropna()
                vals.append(sub.sample(n=min(2500, len(sub)), random_state=10).to_numpy())
                positions.append(pos)
                colors.append(DISPLAY_COLORS[system])
                labels.append(f"{system}\n{rep}")
                pos += 1.0
            pos += 0.6
        bp = ax.boxplot(vals, positions=positions, patch_artist=True, showfliers=False, widths=0.55)
        for patch, color in zip(bp["boxes"], colors):
            patch.set_facecolor(color)
            patch.set_edgecolor("black")
            patch.set_alpha(0.85)
        for median in bp["medians"]:
            median.set_color("black")
            median.set_linewidth(1.1)
        ax.set_xticks(positions)
        ax.set_xticklabels(labels, fontsize=6)
        ax.set_ylabel(ylabel)
        style_axis(ax)
        save_png(fig, OUT_DIRS["microswitch"] / f"replicate_{filename}")


def system_corr(data: pd.DataFrame, system: str, metrics: list[str]) -> pd.DataFrame:
    blocks = []
    for rep in ["rep1", "rep2", "rep3"]:
        cols = {}
        n = None
        for metric in metrics:
            sub = data[(data["system"] == system) & (data["replicate"] == rep) & (data["metric"] == metric)].sort_values("frame")
            vals = sub["value"].to_numpy()
            n = len(vals) if n is None else min(n, len(vals))
            cols[metric] = vals
        if n and n > 0:
            blocks.append(pd.DataFrame({k: v[:n] for k, v in cols.items()}))
    return pd.concat(blocks, ignore_index=True).corr()


def plot_microswitch_correlations() -> None:
    data = pd.read_csv("article/figure/data/microswitch_timeseries.csv")
    metrics = [x[0] for x in MICRO_METRICS]
    labels = [x[3] for x in MICRO_METRICS]
    for system in ["BM213", "C5apep"]:
        fig, ax = plt.subplots(figsize=(4.6, 4.1))
        corr = system_corr(data, system, metrics)
        im = ax.imshow(corr.values, vmin=-1, vmax=1, cmap="RdBu_r")
        ax.set_xticks(range(len(metrics)))
        ax.set_yticks(range(len(metrics)))
        ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=7)
        ax.set_yticklabels(labels, fontsize=7)
        cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
        cbar.set_label("Pearson r", fontsize=7)
        save_png(fig, OUT_DIRS["microswitch"] / f"{system.lower()}_microswitch_correlation.png")


def plot_whole_sequence_rmsf(rmsf_summary: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    for system in SYSTEM_ORDER:
        sub = rmsf_summary[rmsf_summary["system"] == system].sort_values("c5ar_resid")
        x = sub["c5ar_resid"].to_numpy()
        y = sub["mean_rmsf_a"].to_numpy()
        sd = sub["sd_rmsf_a"].fillna(0).to_numpy()
        ax.plot(x, smooth(y, 5), color=DISPLAY_COLORS[system], lw=1.8, label=DISPLAY[system])
        ax.fill_between(x, np.maximum(smooth(y - sd, 5), 0), smooth(y + sd, 5), color=DISPLAY_COLORS[system], alpha=0.16, linewidth=0)
    ax.set_xlabel("C5aR1 residue number")
    ax.set_ylabel("RMSF (Å)")
    ax.set_xlim(20, 325)
    ax.set_ylim(bottom=0)
    ax.legend(fontsize=9.5, loc="upper right")
    style_axis(ax)
    save_png(fig, OUT_DIRS["flexibility"] / "whole_sequence_rmsf_mean_sd.png")


def draw_rmsf_region_track(ax: plt.Axes) -> None:
    ax.set_xlim(20, 325)
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    for label, start, end, color, kind in TM_REGION_BANDS:
        width = end - start
        ax.add_patch(Rectangle((start, 0.28), width, 0.34, facecolor=color, edgecolor="black", linewidth=0.7))
        if kind == "tm":
            ax.text(start + width / 2, 0.75, label, ha="center", va="center", fontsize=10.8, fontweight="bold")
        else:
            ax.text(start + width / 2, 0.13, label, ha="center", va="center", fontsize=10.2, fontweight="bold")


def plot_whole_sequence_rmsf_regions_from_summary(rmsf_summary: pd.DataFrame, output_path: Path) -> None:
    fig = plt.figure(figsize=(8.4, 5.8), constrained_layout=True)
    gs = fig.add_gridspec(2, 1, height_ratios=[3.25, 0.7], hspace=0.05)
    ax = fig.add_subplot(gs[0, 0])
    ax_track = fig.add_subplot(gs[1, 0], sharex=ax)
    for system in SYSTEM_ORDER:
        sub = rmsf_summary[rmsf_summary["system"] == system].sort_values("c5ar_resid")
        x = sub["c5ar_resid"].to_numpy()
        y = sub["mean_rmsf_a"].to_numpy()
        sd = sub["sd_rmsf_a"].fillna(0).to_numpy()
        ax.plot(x, smooth(y, 5), color=DISPLAY_COLORS[system], lw=1.8, label=DISPLAY[system])
        ax.fill_between(x, np.maximum(smooth(y - sd, 5), 0), smooth(y + sd, 5), color=DISPLAY_COLORS[system], alpha=0.14, linewidth=0)
    ax.set_xlabel("C5aR1 residue number")
    ax.set_ylabel("RMSF (Å)")
    ax.set_xlim(20, 325)
    ax.set_ylim(bottom=0)
    ax.legend(fontsize=9.5, loc="upper right")
    style_axis(ax)
    draw_rmsf_region_track(ax_track)
    fig.text(0.02, 0.13, "TM1-7", rotation=90, ha="center", va="center", fontsize=14, fontweight="bold")
    save_png(fig, output_path)


def plot_whole_sequence_rmsf_regions(rmsf_summary: pd.DataFrame) -> None:
    plot_whole_sequence_rmsf_regions_from_summary(rmsf_summary, OUT_DIRS["flexibility"] / "whole_sequence_rmsf_regions.png")


def plot_whole_sequence_rmsf_regions_bm213_no_rep3(rmsf_df: pd.DataFrame) -> None:
    filtered = rmsf_df[~((rmsf_df["system"] == "BM213") & (rmsf_df["replicate"] == "rep3"))].copy()
    rmsf_summary = (
        filtered.groupby(["system", "c5ar_resid"], as_index=False)
        .agg(mean_rmsf_a=("rmsf_a", "mean"), sd_rmsf_a=("rmsf_a", "std"), n_reps=("replicate", "nunique"))
        .sort_values(["system", "c5ar_resid"])
    )
    plot_whole_sequence_rmsf_regions_from_summary(
        rmsf_summary,
        OUT_DIRS["flexibility"] / "whole_sequence_rmsf_regions_bm213_no_rep3.png",
    )


def main() -> None:
    ensure_output_dirs()
    apply_style()
    rmsd_df, rmsf_df, rmsf_summary = build_ca_dynamics()
    plot_rmsd_timeseries(rmsd_df)
    plot_rmsd_distribution(rmsd_df)
    plot_hbond_basic_panels()
    bsa = build_iwi_bsa_replicates()
    plot_iwi_bsa_replicates(bsa)
    plot_iwi_bsa_system_mean(bsa)
    plot_microswitch_by_system()
    plot_microswitch_replicates()
    plot_microswitch_correlations()
    plot_whole_sequence_rmsf(rmsf_summary)
    plot_whole_sequence_rmsf_regions(rmsf_summary)
    plot_whole_sequence_rmsf_regions_bm213_no_rep3(rmsf_df)


if __name__ == "__main__":
    main()
