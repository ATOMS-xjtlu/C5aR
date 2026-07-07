from __future__ import annotations

import pandas as pd
import matplotlib.pyplot as plt

from figure_utils import ROOT, SUPP_DIR, SYSTEM_COLORS, apply_style, image_panel, ensure_dirs, panel_label, save_png, style_axis


def main() -> None:
    ensure_dirs()
    apply_style()
    ci = pd.read_csv("article/figure/data/mdpath_path_confidence.csv")
    nmi = []
    for system, path in [
        ("BM213", ROOT / "gpcr/bm213/mdpath/mdpath_out/nmi_df.csv"),
        ("C5apep", ROOT / "gpcr/c5apep/mdpath/mdpath_out/nmi_df.csv"),
    ]:
        df = pd.read_csv(path).head(15).copy()
        df["system"] = system
        nmi.append(df)
    nmi = pd.concat(nmi, ignore_index=True)

    fig, axes = plt.subplots(3, 2, figsize=(12, 11), constrained_layout=True)
    panels = [
        (axes[0, 0], ROOT / "gpcr/bm213/mdpath/mdpath_out/graph.png", "BM213 MDPath graph", "a"),
        (axes[0, 1], ROOT / "gpcr/c5apep/mdpath/mdpath_out/graph.png", "C5apep MDPath graph", "b"),
        (axes[1, 0], ROOT / "gpcr/bm213/mdpath/mdpath_out/clustered_paths.png", "BM213 clustered paths", "c"),
        (axes[1, 1], ROOT / "gpcr/c5apep/mdpath/mdpath_out/clustered_paths.png", "C5apep clustered paths", "d"),
    ]
    for ax, path, title, label in panels:
        panel_label(ax, label)
        image_panel(ax, path, title)

    ax = axes[2, 0]
    panel_label(ax, "e")
    for system, color in SYSTEM_COLORS.items():
        sub = ci[ci["system"] == system]
        ax.hist(sub["mean"], bins=30, alpha=0.55, color=color, edgecolor="black", label=system)
    ax.set_xlabel("Bootstrap mean path support")
    ax.set_ylabel("Path count")
    ax.legend()
    ax.set_title("Path confidence interval diagnostic")
    style_axis(ax)

    ax = axes[2, 1]
    panel_label(ax, "f")
    for system, color in SYSTEM_COLORS.items():
        sub = nmi[nmi["system"] == system].reset_index(drop=True)
        ax.plot(range(1, len(sub) + 1), sub["MI Difference"], marker="o", color=color, label=system)
    ax.set_xlabel("Top listed NMI pair rank")
    ax.set_ylabel("MI difference")
    ax.legend()
    ax.set_title("Top NMI differences")
    style_axis(ax)

    fig.suptitle("MDPath diagnostics", fontsize=13, fontweight="bold")
    save_png(fig, SUPP_DIR / "FigS5_mdpath_diagnostics.png")


if __name__ == "__main__":
    main()
