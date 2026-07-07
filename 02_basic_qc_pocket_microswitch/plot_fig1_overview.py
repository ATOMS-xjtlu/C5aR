from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from figure_utils import BASIC_MAIN_DIR, MAIN_DIR, SYSTEM_COLORS, apply_style, draw_placeholder, ensure_dirs, grouped_bar_points, panel_label, save_png, save_pngs, style_axis


def plot_contact_panel(ax, summary: pd.DataFrame, metric: str, ylabel: str, ylim: tuple[float, float], title: str, label: str) -> None:
    panel_label(ax, label)
    grouped_bar_points(ax, summary, metric, ylabel, ylim)
    ax.set_title(title, fontsize=10)


def save_contact_panel(summary: pd.DataFrame, metric: str, ylabel: str, ylim: tuple[float, float], title: str, label: str, output_path) -> None:
    fig, ax = plt.subplots(figsize=(4.4, 4.8))
    plot_contact_panel(ax, summary, metric, ylabel, ylim, title, label)
    save_png(fig, output_path)


def main() -> None:
    ensure_dirs()
    apply_style()
    hbond = pd.read_csv("article/figure/tables/hbond_summary.csv")
    basin = pd.read_csv("article/figure/tables/fes_basin_summary.csv")
    mdpath = pd.read_csv("article/figure/tables/mdpath_summary.csv")

    fig = plt.figure(figsize=(11, 7.2))
    gs = fig.add_gridspec(2, 3, width_ratios=[1.35, 1, 1], height_ratios=[1.0, 1.0], wspace=0.38, hspace=0.42)

    ax = fig.add_subplot(gs[:, 0])
    panel_label(ax, "a")
    ax.axis("off")
    boxes = [
        ("Orthosteric pocket", "Y258 vs S171 contacts\nIWI pocket engagement"),
        ("Core microswitches", "M120/W255/R3.50\nTM6-TM7 packing"),
        ("Allosteric network", "FES projection\nMDPath communication"),
    ]
    y_positions = [0.78, 0.50, 0.22]
    for (title, body), y in zip(boxes, y_positions):
        rect = plt.Rectangle((0.08, y - 0.09), 0.82, 0.16, facecolor="#F7F7F7", edgecolor="black", lw=0.9, transform=ax.transAxes)
        ax.add_patch(rect)
        ax.text(0.12, y + 0.025, title, transform=ax.transAxes, fontweight="bold", fontsize=10, va="center")
        ax.text(0.12, y - 0.045, body, transform=ax.transAxes, fontsize=8, va="center")
    for y1, y2 in zip(y_positions[:-1], y_positions[1:]):
        ax.annotate("", xy=(0.49, y2 + 0.09), xytext=(0.49, y1 - 0.09), xycoords=ax.transAxes, arrowprops=dict(arrowstyle="-|>", lw=1.1, color="#444444"))
    ax.text(0.08, 0.04, "Main scope: BM213 vs C5apep, three replicates each\nOutput format: PNG only", transform=ax.transAxes, fontsize=8)

    ax = fig.add_subplot(gs[0, 1])
    plot_contact_panel(ax, hbond, "Y258_hbond", "Mean H-bond count", (0, 1.1), "BM213 stabilizes Y258 contact", "b")

    ax = fig.add_subplot(gs[0, 2])
    plot_contact_panel(ax, hbond, "S171_hbond", "Mean H-bond count", (0, 1.8), "C5apep favors S171 contact", "c")

    ax = fig.add_subplot(gs[1, 1])
    panel_label(ax, "d")
    occ = basin.groupby(["system", "basin"], as_index=False)["occupancy"].mean()
    pivot = occ.pivot(index="basin", columns="system", values="occupancy").fillna(0)
    x = np.arange(len(pivot))
    width = 0.34
    for i, system in enumerate(["BM213", "C5apep"]):
        ax.bar(x + (i - 0.5) * width, pivot.get(system, 0), width=width, color=SYSTEM_COLORS[system], edgecolor="black", label=system)
    ax.set_xticks(x)
    ax.set_xticklabels(pivot.index, rotation=20, ha="right")
    ax.set_ylabel("Mean occupancy")
    ax.set_ylim(0, 1)
    ax.legend(loc="upper right")
    ax.set_title("FES projection occupancy", fontsize=10)
    style_axis(ax)

    ax = fig.add_subplot(gs[1, 2])
    panel_label(ax, "e")
    weights = mdpath.set_index("system")["top_weight"]
    ax.bar(weights.index, weights.values, color=[SYSTEM_COLORS[s] for s in weights.index], edgecolor="black", width=0.6)
    ax.set_ylabel("Top path weight")
    ax.set_ylim(0, max(weights.values) * 1.25)
    ax.set_title("MDPath top communication path", fontsize=10)
    style_axis(ax)

    fig.suptitle("C5aR1 ligand-dependent signaling model", y=0.995, fontsize=13, fontweight="bold")
    save_pngs(fig, [MAIN_DIR / "Fig1_overview.png", BASIC_MAIN_DIR / "Fig1_overview.png"])
    save_contact_panel(hbond, "Y258_hbond", "Mean H-bond count", (0, 1.1), "BM213 stabilizes Y258 contact", "b", BASIC_MAIN_DIR / "Fig1_overview_b.png")
    save_contact_panel(hbond, "S171_hbond", "Mean H-bond count", (0, 1.8), "C5apep favors S171 contact", "c", BASIC_MAIN_DIR / "Fig1_overview_c.png")


if __name__ == "__main__":
    main()
