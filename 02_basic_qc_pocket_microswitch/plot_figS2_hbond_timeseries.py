from __future__ import annotations

import pandas as pd
import matplotlib.pyplot as plt

from figure_utils import SUPP_DIR, SYSTEM_COLORS, apply_style, downsample, ensure_dirs, panel_label, rolling_mean, save_png, style_axis


def main() -> None:
    ensure_dirs()
    apply_style()
    ts = pd.read_csv("article/figure/data/hbond_timeseries.csv")
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.2), sharex=True, constrained_layout=True)
    panels = [
        ("BM213", "Y258_hbond", "a"),
        ("C5apep", "Y258_hbond", "b"),
        ("BM213", "S171_hbond", "c"),
        ("C5apep", "S171_hbond", "d"),
    ]
    for ax, (system, metric, label) in zip(axes.ravel(), panels):
        panel_label(ax, label)
        sub0 = ts[(ts["system"] == system) & (ts["metric"] == metric)]
        for rep, sub in sub0.groupby("replicate"):
            sub = sub.sort_values("frame").copy()
            sub["smooth"] = rolling_mean(sub["value"], 180)
            plot = downsample(sub, 1800)
            ax.plot(plot["time_ns"], plot["smooth"], lw=1.1, alpha=0.85, label=rep)
        ax.set_title(f"{system} {metric.replace('_', ' ')}")
        ax.set_ylabel("Rolling mean count")
        ax.set_xlim(0, 2000)
        ax.legend(fontsize=7)
        style_axis(ax)
    axes[1, 0].set_xlabel("Time (ns)")
    axes[1, 1].set_xlabel("Time (ns)")
    fig.suptitle("Y258 and S171 H-bond time series across replicates", fontsize=13, fontweight="bold")
    save_png(fig, SUPP_DIR / "FigS2_hbond_timeseries.png")


if __name__ == "__main__":
    main()
