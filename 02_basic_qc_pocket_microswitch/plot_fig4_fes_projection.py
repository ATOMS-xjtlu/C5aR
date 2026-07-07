from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from figure_utils import MAIN_DIR, SYSTEM_COLORS, apply_style, draw_placeholder, ensure_dirs, make_free_energy, panel_label, save_png, style_axis


def fes_panel(ax, data: pd.DataFrame, system: str, ranges):
    sub = data[data["system"] == system]
    fes, xedges, yedges = make_free_energy(sub["M120_distance"].to_numpy(), sub["DRY_PIF_distance"].to_numpy(), bins=60, ranges=ranges)
    cmap = plt.get_cmap("viridis_r").copy()
    cmap.set_bad("#ECECEC")
    im = ax.imshow(fes, origin="lower", aspect="auto", extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], cmap=cmap, vmin=0, vmax=np.nanpercentile(fes, 90))
    ax.contour((xedges[:-1] + xedges[1:]) / 2, (yedges[:-1] + yedges[1:]) / 2, fes, levels=[1, 2, 3, 4], colors="white", linewidths=0.45, alpha=0.8)
    ax.set_title(system)
    ax.set_xlabel("M120-D82 distance (Å)")
    ax.set_ylabel("Y222$^{5.58}$-F251$^{6.44}$ distance (Å)")
    style_axis(ax)
    return im


def main() -> None:
    ensure_dirs()
    apply_style()
    data = pd.read_csv("article/figure/data/fes_points.csv")
    basin = pd.read_csv("article/figure/tables/fes_basin_summary.csv")
    xlim = (data["M120_distance"].quantile(0.005), data["M120_distance"].quantile(0.995))
    ylim = (data["DRY_PIF_distance"].quantile(0.005), data["DRY_PIF_distance"].quantile(0.995))

    fig, axes = plt.subplots(2, 2, figsize=(10.8, 8.0), constrained_layout=True)
    ax = axes[0, 0]
    panel_label(ax, "a")
    im = fes_panel(ax, data, "BM213", [xlim, ylim])
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    cbar.set_label("Free energy (kcal/mol)", fontsize=8)

    ax = axes[0, 1]
    panel_label(ax, "b")
    im = fes_panel(ax, data, "C5apep", [xlim, ylim])
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    cbar.set_label("Free energy (kcal/mol)", fontsize=8)

    ax = axes[1, 0]
    panel_label(ax, "c")
    occ = basin.groupby(["system", "basin"], as_index=False)["occupancy"].mean()
    pivot = occ.pivot(index="basin", columns="system", values="occupancy").fillna(0)
    x = np.arange(len(pivot))
    width = 0.34
    for i, system in enumerate(["BM213", "C5apep"]):
        ax.bar(x + (i - 0.5) * width, pivot.get(system, 0), width=width, color=SYSTEM_COLORS[system], edgecolor="black", label=system)
    ax.set_xticks(x)
    ax.set_xticklabels(pivot.index)
    ax.set_ylabel("Mean occupancy")
    ax.set_ylim(0, 1)
    ax.legend()
    ax.set_title("Data-driven projection basins")
    style_axis(ax)

    ax = axes[1, 1]
    panel_label(ax, "d")
    draw_placeholder(
        ax,
        "Representative basin structure placeholder",
        [
            "Use article/structure/02_representative_states/",
            "Recommended: BM213/C5apep TM67 c0-c2",
            "Final replacement should be a PyMOL JPG/PNG render",
        ],
    )

    fig.suptitle("Free-energy projection of microswitch states", fontsize=13, fontweight="bold")
    save_png(fig, MAIN_DIR / "Fig4_fes_projection.png")


if __name__ == "__main__":
    main()
