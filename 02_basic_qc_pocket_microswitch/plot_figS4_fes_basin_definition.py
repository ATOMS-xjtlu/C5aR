from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from figure_utils import SUPP_DIR, apply_style, downsample, ensure_dirs, panel_label, save_png, style_axis

BASIN_COLORS = {"Basin 1": "#66C2A5", "Basin 2": "#8DA0CB", "Basin 3": "#FC8D62"}


def main() -> None:
    ensure_dirs()
    apply_style()
    data = pd.read_csv("article/figure/data/fes_points.csv")
    basin = pd.read_csv("article/figure/tables/fes_basin_summary.csv")
    fig, axes = plt.subplots(2, 2, figsize=(10.8, 7.8), constrained_layout=True)

    for ax, system, label in [(axes[0, 0], "BM213", "a"), (axes[0, 1], "C5apep", "b")]:
        panel_label(ax, label)
        sub = downsample(data[data["system"] == system], 7000)
        for basin_name, bsub in sub.groupby("basin"):
            ax.scatter(bsub["M120_distance"], bsub["DRY_PIF_distance"], s=2, alpha=0.28, color=BASIN_COLORS.get(basin_name, "#999999"), label=basin_name)
        ax.set_xlabel("M120-D82 distance (Å)")
        ax.set_ylabel("Y222$^{5.58}$-F251$^{6.44}$ distance (Å)")
        ax.set_title(f"{system} projection basin assignment")
        ax.legend(markerscale=4, fontsize=7)
        style_axis(ax)

    ax = axes[1, 0]
    panel_label(ax, "c")
    mean_occ = basin.groupby(["system", "basin"], as_index=False)["occupancy"].mean()
    for i, system in enumerate(["BM213", "C5apep"]):
        sub = mean_occ[mean_occ["system"] == system].sort_values("basin")
        ax.bar(np.arange(len(sub)) + i * 0.35, sub["occupancy"], width=0.34, label=system, edgecolor="black", color=["#66C2A5", "#8DA0CB", "#FC8D62"][: len(sub)], alpha=0.85 if system == "BM213" else 0.55)
    ax.set_xticks(np.arange(3) + 0.175)
    ax.set_xticklabels(["Basin 1", "Basin 2", "Basin 3"])
    ax.set_ylabel("Mean occupancy")
    ax.set_ylim(0, 1)
    ax.legend()
    ax.set_title("Basin occupancy")
    style_axis(ax)

    ax = axes[1, 1]
    panel_label(ax, "d")
    ax.axis("off")
    centroids = basin.groupby("basin", as_index=False)[["centroid_M120_distance", "centroid_DRY_PIF_distance"]].mean()
    lines = ["Basin    M120-D82 centroid   Y222-F251 centroid"]
    for _, row in centroids.iterrows():
        lines.append(f"{row.basin:8s} {row.centroid_M120_distance:13.2f} {row.centroid_DRY_PIF_distance:17.2f}")
    lines.extend(["", "Definition: three data-driven clusters", "in the M120-D82 vs Y222-F251 projection."])
    ax.text(0.05, 0.92, "\n".join(lines), transform=ax.transAxes, va="top", family="monospace", fontsize=8)
    ax.set_title("Basin definition summary")

    fig.suptitle("FES/projection basin definition and occupancy", fontsize=13, fontweight="bold")
    save_png(fig, SUPP_DIR / "FigS4_fes_basin_definition.png")


if __name__ == "__main__":
    main()
