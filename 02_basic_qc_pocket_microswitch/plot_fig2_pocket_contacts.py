from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from figure_utils import BASIC_MAIN_DIR, DATA_DIR, MAIN_DIR, SYSTEM_COLORS, TABLE_DIR, apply_style, downsample, draw_placeholder, ensure_dirs, grouped_bar_points, panel_label, rolling_mean, save_png, save_pngs, style_axis


IWI_BSA_INPUTS = {
    "BM213": {
        "label": "BM213 rep1",
        "sc": Path("gpcr/bm213/BSA_IWI/sasa_sc_58_69_83.dat"),
        "lig": Path("gpcr/bm213/BSA_IWI/sasa_lig.dat"),
        "total": Path("gpcr/bm213/BSA_IWI/sasa_total.dat"),
    },
    "C5apep": {
        "label": "C5apep (PEP) rep3",
        "sc": Path("gpcr/c5apep/third/BSA_IWI/sasa_sc_60_71_85.dat"),
        "lig": Path("gpcr/c5apep/third/BSA_IWI/sasa_lig_284.dat"),
        "total": Path("gpcr/c5apep/third/BSA_IWI/sasa_total.dat"),
    },
}


def read_sasa(path: Path, column: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep=r"\s+", comment="#", names=["frame", column])
    df["frame"] = df["frame"].astype(int)
    return df


def build_iwi_bsa() -> pd.DataFrame:
    blocks = []
    for system, spec in IWI_BSA_INPUTS.items():
        sc = read_sasa(spec["sc"], "sasa_iwi")
        lig = read_sasa(spec["lig"], "sasa_ligand")
        total = read_sasa(spec["total"], "sasa_complex")
        df = sc.merge(lig, on="frame").merge(total, on="frame")
        df["system"] = system
        df["system_label"] = spec["label"]
        df["buried_sasa_a2"] = (df["sasa_iwi"] + df["sasa_ligand"] - df["sasa_complex"]).clip(lower=0)
        max_frame = max(int(df["frame"].max()), 1)
        df["time_ns"] = (df["frame"] - 1) * (2000.0 / (max_frame - 1)) if max_frame > 1 else 0.0
        blocks.append(df[["system", "system_label", "frame", "time_ns", "sasa_iwi", "sasa_ligand", "sasa_complex", "buried_sasa_a2"]])
    out = pd.concat(blocks, ignore_index=True)
    out.to_csv(DATA_DIR / "iwi_bsa_timeseries.csv", index=False)
    summary = (
        out.groupby(["system", "system_label"], as_index=False)["buried_sasa_a2"]
        .agg(frames="count", mean="mean", sd="std", median="median", q25=lambda x: x.quantile(0.25), q75=lambda x: x.quantile(0.75), minimum="min", maximum="max")
        .round(5)
    )
    summary.to_csv(TABLE_DIR / "iwi_bsa_summary.csv", index=False)
    return out


def plot_iwi_bsa_panel(ax: plt.Axes, bsa: pd.DataFrame, label: str) -> None:
    panel_label(ax, label)
    systems = ["BM213", "C5apep"]
    vals = []
    for system in systems:
        sub = bsa[bsa["system"] == system]["buried_sasa_a2"].dropna()
        vals.append(sub.sample(n=min(5000, len(sub)), random_state=8).to_numpy())
    bp = ax.boxplot(vals, patch_artist=True, widths=0.58, showfliers=False)
    for patch, system in zip(bp["boxes"], systems):
        patch.set_facecolor(SYSTEM_COLORS[system])
        patch.set_alpha(0.85)
        patch.set_edgecolor("black")
    for median in bp["medians"]:
        median.set_color("black")
        median.set_linewidth(1.2)
    means = [float(np.mean(v)) for v in vals]
    ax.scatter([1, 2], means, s=34, marker="D", color="white", edgecolor="black", zorder=4, label="mean")
    ax.set_xticks([1, 2])
    ax.set_xticklabels(["BM213\nrep1", "C5apep\n(PEP rep3)"])
    ax.set_ylabel("IWI buried SASA (Å²)")
    ax.set_title("IWI pocket burial", fontsize=10)
    ax.text(0.03, 0.96, "I91/W102/I116 + ligand", transform=ax.transAxes, va="top", fontsize=7)
    ax.legend(fontsize=7, loc="upper right")
    style_axis(ax)


def save_iwi_bsa_panel(bsa: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(4.6, 4.8))
    plot_iwi_bsa_panel(ax, bsa, "e")
    save_png(fig, BASIC_MAIN_DIR / "Fig2_iwi_bsa.png")


def main() -> None:
    ensure_dirs()
    apply_style()
    ts = pd.read_csv("article/figure/data/hbond_timeseries.csv")
    summary = pd.read_csv("article/figure/tables/hbond_summary.csv")
    iwi_bsa = build_iwi_bsa()

    fig = plt.figure(figsize=(13.2, 7.4), constrained_layout=True)
    gs = fig.add_gridspec(2, 3, width_ratios=[1.08, 1.0, 1.0], height_ratios=[1, 1])
    ax = fig.add_subplot(gs[0, 0])
    panel_label(ax, "a")
    draw_placeholder(
        ax,
        "Pocket close-up placeholder",
        [
            "Replace with PyMOL render from:",
            "article/structure/01_pocket_closeup/",
            "Labels: Y258, S171, W255, M120, ligand",
        ],
    )

    ax = fig.add_subplot(gs[0, 1])
    panel_label(ax, "b")
    grouped_bar_points(ax, summary, "Y258_hbond", "Mean H-bond count", (0, 1.1))
    ax.set_title("Y258-ligand contact")

    ax = fig.add_subplot(gs[0, 2])
    panel_label(ax, "c")
    grouped_bar_points(ax, summary, "S171_hbond", "Mean H-bond count", (0, 1.8))
    ax.set_title("S171-ligand contact")

    ax = fig.add_subplot(gs[1, :2])
    panel_label(ax, "d")
    for system, color in SYSTEM_COLORS.items():
        for metric, ls in [("Y258_hbond", "-"), ("S171_hbond", "--")]:
            sub = ts[(ts["system"] == system) & (ts["metric"] == metric)]
            mean_ts = sub.groupby("frame", as_index=False).agg(time_ns=("time_ns", "mean"), value=("value", "mean"))
            mean_ts["smooth"] = rolling_mean(mean_ts["value"], 300)
            plot = downsample(mean_ts, 1800)
            label = f"{system} {metric.split('_')[0]}"
            ax.plot(plot["time_ns"], plot["smooth"], color=color, lw=1.5, ls=ls, label=label)
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Rolling mean H-bond count")
    ax.set_xlim(0, 2000)
    ax.legend(fontsize=7, ncol=2, loc="upper right")
    style_axis(ax)

    ax = fig.add_subplot(gs[1, 2])
    plot_iwi_bsa_panel(ax, iwi_bsa, "e")

    fig.suptitle("Ligand-specific orthosteric pocket contacts", fontsize=13, fontweight="bold")
    save_pngs(fig, [MAIN_DIR / "Fig2_pocket_contacts.png", BASIC_MAIN_DIR / "Fig2_pocket_contacts.png"])
    save_iwi_bsa_panel(iwi_bsa)


if __name__ == "__main__":
    main()
