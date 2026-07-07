from __future__ import annotations

import pandas as pd
import matplotlib.pyplot as plt

from figure_utils import BASIC_SUPP_DIR, METRIC_LABELS, SUPP_DIR, SYSTEM_COLORS, apply_style, ensure_dirs, panel_label, save_png, save_pngs, style_axis


def plot_metric(ax: plt.Axes, data: pd.DataFrame, metric: str, label: str) -> None:
    panel_label(ax, label)
    vals, positions, colors, xticks = [], [], [], []
    pos = 1
    for system in ["BM213", "C5apep"]:
        for rep in ["rep1", "rep2", "rep3"]:
            sub = data[(data["system"] == system) & (data["replicate"] == rep) & (data["metric"] == metric)]["value"].dropna()
            vals.append(sub.sample(n=min(2500, len(sub)), random_state=10).to_numpy())
            positions.append(pos)
            colors.append(SYSTEM_COLORS[system])
            xticks.append(f"{system}\n{rep}")
            pos += 1
        pos += 0.6
    bp = ax.boxplot(vals, positions=positions, patch_artist=True, showfliers=False, widths=0.55)
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_edgecolor("black")
        patch.set_alpha(0.85)
    ax.set_xticks(positions)
    ax.set_xticklabels(xticks, fontsize=6)
    ax.set_ylabel(METRIC_LABELS.get(metric, metric))
    style_axis(ax)


def save_single_metric(data: pd.DataFrame, metric: str, letter: str) -> None:
    fig, ax = plt.subplots(figsize=(5.0, 4.8))
    plot_metric(ax, data, metric, letter)
    fig.suptitle(f"{METRIC_LABELS.get(metric, metric)} across replicates", fontsize=12, fontweight="bold")
    save_png(fig, BASIC_SUPP_DIR / f"FigS3_microswitch_replicates_{letter}.png")


def main() -> None:
    ensure_dirs()
    apply_style()
    data = pd.read_csv("article/figure/data/microswitch_timeseries.csv")
    metrics = ["M120_distance", "DRY_PIF_distance", "W255_chi2", "R350_chi2", "TM62_IC_distance", "TM72_IC_distance"]
    fig, axes = plt.subplots(2, 3, figsize=(12, 7.4), constrained_layout=True)
    labels = list("abcdef")
    for ax, metric, label in zip(axes.ravel(), metrics, labels):
        plot_metric(ax, data, metric, label)
    fig.suptitle("Microswitch metric distributions for all replicates", fontsize=13, fontweight="bold")
    save_pngs(fig, [SUPP_DIR / "FigS3_microswitch_replicates.png", BASIC_SUPP_DIR / "FigS3_microswitch_replicates.png"])

    for metric, label in zip(metrics, labels):
        save_single_metric(data, metric, label)


if __name__ == "__main__":
    main()
