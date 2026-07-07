from __future__ import annotations

import pandas as pd
import matplotlib.pyplot as plt

from figure_utils import MAIN_DIR, ROOT, SYSTEM_COLORS, apply_style, draw_placeholder, ensure_dirs, image_panel, panel_label, save_png, style_axis


def read_cluster(path):
    return pd.read_csv(path, sep=r"\s+", comment="#", names=["Cluster", "Frames", "Frac", "AvgDist", "Stdev", "Centroid", "AvgCDist"])


def main() -> None:
    ensure_dirs()
    apply_style()
    paths = pd.read_csv("article/figure/data/mdpath_top_paths.csv")
    summary = pd.read_csv("article/figure/tables/mdpath_summary.csv")

    fig, axes = plt.subplots(2, 2, figsize=(11, 7.6), constrained_layout=True)

    ax = axes[0, 0]
    panel_label(ax, "a")
    for system, color in SYSTEM_COLORS.items():
        sub = paths[paths["system"] == system].head(15)
        offset = -0.18 if system == "BM213" else 0.18
        ax.bar(sub["rank"] + offset, sub["total_weight"], width=0.34, color=color, edgecolor="black", label=system)
    ax.set_xlabel("Top path rank")
    ax.set_ylabel("Total weight")
    ax.legend()
    ax.set_title("MDPath top communication paths")
    style_axis(ax)

    ax = axes[0, 1]
    panel_label(ax, "b")
    x = range(len(summary))
    ax.bar([i - 0.16 for i in x], summary["top_weight"], width=0.32, color=[SYSTEM_COLORS[s] for s in summary["system"]], edgecolor="black", label="Top weight")
    ax.bar([i + 0.16 for i in x], summary["mean_top20_weight"], width=0.32, color="#BBBBBB", edgecolor="black", label="Mean top 20")
    ax.set_xticks(list(x))
    ax.set_xticklabels(summary["system"])
    ax.set_ylabel("Path weight")
    ax.legend(fontsize=7)
    ax.set_title("Path strength summary")
    style_axis(ax)

    ax = axes[1, 0]
    panel_label(ax, "c")
    bm = read_cluster(ROOT / "article/structure/04_reference_mapping/BM213_TM67_cluster_summary.dat").head(5)
    pep = read_cluster(ROOT / "article/structure/04_reference_mapping/C5apep_TM67_cluster_summary.dat").head(5)
    ax.plot(bm["Cluster"], bm["Frac"], marker="o", lw=1.8, color=SYSTEM_COLORS["BM213"], label="BM213 TM67")
    ax.plot(pep["Cluster"], pep["Frac"], marker="o", lw=1.8, color=SYSTEM_COLORS["C5apep"], label="C5apep TM67")
    ax.set_xlabel("Cluster")
    ax.set_ylabel("Fraction")
    ax.set_ylim(0, max(bm["Frac"].max(), pep["Frac"].max()) * 1.25)
    ax.legend()
    ax.set_title("Representative TM6/TM7 cluster fractions")
    style_axis(ax)

    ax = axes[1, 1]
    panel_label(ax, "d")
    draw_placeholder(
        ax,
        "MDPath 3D pathway placeholder",
        [
            "Use article/structure/03_mdpath_3d/",
            "BM213/C5apep quick or full PyMOL loaders",
            "Replace this panel with manually rendered 3D pathway view",
        ],
    )

    fig.suptitle("MDPath network and representative pathway summary", fontsize=13, fontweight="bold")
    save_png(fig, MAIN_DIR / "Fig5_mdpath_network.png")


if __name__ == "__main__":
    main()
