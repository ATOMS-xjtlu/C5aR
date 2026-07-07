from __future__ import annotations

import pandas as pd
import matplotlib.pyplot as plt

from figure_utils import BASIC_MAIN_DIR, MAIN_DIR, SYSTEM_COLORS, apply_style, ensure_dirs, panel_label, plot_kde_or_hist, save_png, save_pngs, style_axis


PANEL_A_LABEL = "M120-D82 distance (Å)"
PANEL_B_LABEL = "W255 χ2 (deg)"
PANEL_C_LABEL = "R3.50 χ2 (deg)"
PANEL_D_LABEL = "TM6/TM2 and TM7/TM2 distance (Å)"
SHORT_LABELS = ["M120-D82", "Y222-F251", "W255", "R3.50", "TM6/TM2", "TM7/TM2"]


def box_by_system(ax: plt.Axes, data: pd.DataFrame, metric: str, ylabel: str) -> None:
    vals = [
        data[(data["system"] == s) & (data["metric"] == metric)]["value"]
        .dropna()
        .sample(n=min(5000, len(data[(data["system"] == s) & (data["metric"] == metric)])), random_state=4)
        .to_numpy()
        for s in ["BM213", "C5apep"]
    ]
    bp = ax.boxplot(vals, patch_artist=True, widths=0.55, showfliers=False)
    for patch, system in zip(bp["boxes"], ["BM213", "C5apep"]):
        patch.set_facecolor(SYSTEM_COLORS[system])
        patch.set_alpha(0.85)
        patch.set_edgecolor("black")
    for item in bp["medians"]:
        item.set_color("black")
        item.set_linewidth(1.2)
    ax.set_xticks([1, 2])
    ax.set_xticklabels(["BM213", "C5apep"])
    ax.set_ylabel(ylabel)
    style_axis(ax)


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


def plot_m120_panel(ax: plt.Axes, data: pd.DataFrame, label: str) -> None:
    panel_label(ax, label)
    box_by_system(ax, data, "M120_distance", PANEL_A_LABEL)
    ax.set_title("M120 vs D82", fontsize=10)


def plot_chi_panel(ax: plt.Axes, data: pd.DataFrame, metric: str, label: str, title: str, xlabel: str) -> None:
    panel_label(ax, label)
    for system, color in SYSTEM_COLORS.items():
        sub = data[(data["system"] == system) & (data["metric"] == metric)]
        plot_kde_or_hist(ax, sub["value"].to_numpy(), color, system, bins=90)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Density")
    ax.legend()
    ax.set_title(title, fontsize=10)
    style_axis(ax)


def plot_tm_box_panel(ax: plt.Axes, data: pd.DataFrame, label: str) -> None:
    panel_label(ax, label)
    combo = data[data["metric"].isin(["TM62_IC_distance", "TM72_IC_distance"])].copy()
    positions, vals, colors, labels = [], [], [], []
    pos = 1
    for metric, metric_label in [("TM62_IC_distance", "TM6/TM2"), ("TM72_IC_distance", "TM7/TM2")]:
        for system in ["BM213", "C5apep"]:
            sub = combo[(combo["metric"] == metric) & (combo["system"] == system)]["value"].dropna()
            vals.append(sub.sample(n=min(5000, len(sub)), random_state=5).to_numpy())
            positions.append(pos)
            colors.append(SYSTEM_COLORS[system])
            labels.append(f"{metric_label}\n{system}")
            pos += 1
        pos += 0.45
    bp = ax.boxplot(vals, positions=positions, patch_artist=True, showfliers=False, widths=0.5)
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_edgecolor("black")
    ax.set_xticks(positions)
    ax.set_xticklabels(labels, fontsize=7)
    ax.set_ylabel("Distance (Å)")
    ax.set_title("TM6/TM7 intracellular opening", fontsize=10)
    style_axis(ax)


def plot_corr_panel(ax: plt.Axes, data: pd.DataFrame, system: str, label: str) -> None:
    panel_label(ax, label)
    metrics = ["M120_distance", "DRY_PIF_distance", "W255_chi2", "R350_chi2", "TM62_IC_distance", "TM72_IC_distance"]
    corr = system_corr(data, system, metrics)
    im = ax.imshow(corr.values, vmin=-1, vmax=1, cmap="RdBu_r")
    ax.set_xticks(range(len(metrics)))
    ax.set_yticks(range(len(metrics)))
    ax.set_xticklabels(SHORT_LABELS, rotation=45, ha="right", fontsize=7)
    ax.set_yticklabels(SHORT_LABELS, fontsize=7)
    ax.set_title(f"{system} metric correlation", fontsize=10)
    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    cbar.set_label("Pearson r", fontsize=7)


def save_single_panel(data: pd.DataFrame, output_path, render_fn, figsize=(4.4, 4.8)) -> None:
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)
    render_fn(ax)
    save_png(fig, output_path)


def main() -> None:
    ensure_dirs()
    apply_style()
    data = pd.read_csv("article/figure/data/microswitch_timeseries.csv")
    fig, axes = plt.subplots(2, 3, figsize=(12, 7.5), constrained_layout=True)

    plot_m120_panel(axes[0, 0], data, "a")
    plot_chi_panel(axes[0, 1], data, "W255_chi2", "b", "W255 side-chain rotamer", PANEL_B_LABEL)
    plot_chi_panel(axes[0, 2], data, "R350_chi2", "c", "R3.50 side-chain rotamer", PANEL_C_LABEL)
    plot_tm_box_panel(axes[1, 0], data, "d")
    plot_corr_panel(axes[1, 1], data, "BM213", "e")
    plot_corr_panel(axes[1, 2], data, "C5apep", "f")

    fig.suptitle("Core microswitch and intracellular packing rearrangements", fontsize=13, fontweight="bold")
    save_pngs(fig, [MAIN_DIR / "Fig3_microswitch.png", BASIC_MAIN_DIR / "Fig3_microswitch.png"])

    single_panels = [
        ("a", lambda ax: plot_m120_panel(ax, data, "a"), (4.4, 4.8)),
        ("b", lambda ax: plot_chi_panel(ax, data, "W255_chi2", "b", "W255 side-chain rotamer", PANEL_B_LABEL), (4.6, 4.8)),
        ("c", lambda ax: plot_chi_panel(ax, data, "R350_chi2", "c", "R3.50 side-chain rotamer", PANEL_C_LABEL), (4.6, 4.8)),
        ("d", lambda ax: plot_tm_box_panel(ax, data, "d"), (5.2, 4.8)),
        ("e", lambda ax: plot_corr_panel(ax, data, "BM213", "e"), (5.0, 4.8)),
        ("f", lambda ax: plot_corr_panel(ax, data, "C5apep", "f"), (5.0, 4.8)),
    ]
    for letter, render_fn, figsize in single_panels:
        save_single_panel(data, BASIC_MAIN_DIR / f"Fig3_microswitch_{letter}.png", render_fn, figsize=figsize)


if __name__ == "__main__":
    main()
