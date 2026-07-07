from __future__ import annotations

from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt

from figure_utils import BASIC_SUPP_DIR, SUPP_DIR, SYSTEM_COLORS, add_time_ns, apply_style, downsample, ensure_dirs, panel_label, read_dat, rolling_mean, save_pngs, style_axis


def main() -> None:
    ensure_dirs()
    apply_style()
    rmsd = pd.read_csv("article/figure/data/rmsd_timeseries.csv")
    qc = pd.read_csv("article/figure/tables/production_qc.csv")
    apo_rmsd = add_time_ns(read_dat(Path("gpcr/apo/rmsd_TM56.dat")))

    fig, axes = plt.subplots(2, 2, figsize=(11, 7.5), constrained_layout=True)

    ax = axes[0, 0]
    panel_label(ax, "a")
    for (system, rep), sub in rmsd.groupby(["system", "replicate"]):
        sub = sub.sort_values("frame").copy()
        sub["smooth"] = rolling_mean(sub["value"], 250)
        plot = downsample(sub, 1600)
        ax.plot(plot["time_ns"], plot["smooth"], color=SYSTEM_COLORS[system], lw=1.1, alpha=0.55, label=f"{system} {rep}")
    apo_plot = apo_rmsd.copy()
    apo_plot["smooth"] = rolling_mean(apo_plot["value"], 250)
    apo_plot = downsample(apo_plot, 1600)
    ax.plot(apo_plot["time_ns"], apo_plot["smooth"], color="#444444", lw=1.5, ls="--", label="apo TM5/6 core")
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("RMSD (Å)")
    ax.set_xlim(0, 2000)
    ax.legend(fontsize=6, ncol=3, loc="upper left")
    ax.set_title("Ligand and apo RMSD traces")
    style_axis(ax)

    ax = axes[0, 1]
    panel_label(ax, "b")
    for i, system in enumerate(["BM213", "C5apep"]):
        sub = qc[qc["system"] == system]
        x = range(i * 8, i * 8 + len(sub))
        ax.scatter(list(x), sub["mean_temp_k"], color=SYSTEM_COLORS[system], edgecolor="black", s=30, label=system)
    ax.axhline(303.15, color="#444444", lw=1, ls="--")
    ax.set_ylabel("Mean sampled temperature (K)")
    ax.set_xticks([])
    ax.legend()
    ax.set_title("Production temperature QC")
    style_axis(ax)

    ax = axes[1, 0]
    panel_label(ax, "c")
    for i, system in enumerate(["BM213", "C5apep"]):
        sub = qc[qc["system"] == system]
        x = range(i * 8, i * 8 + len(sub))
        ax.scatter(list(x), sub["mean_density"], color=SYSTEM_COLORS[system], edgecolor="black", s=30, label=system)
    ax.set_ylabel("Mean sampled density")
    ax.set_xticks([])
    ax.legend()
    ax.set_title("Production density QC")
    style_axis(ax)

    ax = axes[1, 1]
    panel_label(ax, "d")
    ax.axis("off")
    counts = qc.groupby(["system", "replicate"]).agg(n_logs=("log_path", "count"), warnings=("error_or_warning_count", "sum"), last_time_ps=("last_time_ps", "max")).reset_index()
    lines = ["System   Rep   logs   warnings   last time (ns)"]
    for _, row in counts.iterrows():
        lines.append(f"{row.system:7s} {row.replicate:4s} {int(row.n_logs):4d} {int(row.warnings):10d} {row.last_time_ps/1000:12.1f}")
    ax.text(0.02, 0.95, "\n".join(lines), transform=ax.transAxes, va="top", family="monospace", fontsize=8)
    ax.set_title("Parsed production-log inventory")

    fig.suptitle("Trajectory and production-log QC", fontsize=13, fontweight="bold")
    save_pngs(fig, [SUPP_DIR / "FigS1_qc.png", BASIC_SUPP_DIR / "FigS1_qc.png"])


if __name__ == "__main__":
    main()
