from __future__ import annotations

from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

from figure_utils import BASIC_SUPP_DIR, apply_style, ensure_dirs, save_png, style_axis


PALETTE = {
    "apo": "#6B7280",
    "c5a": "#D96B6B",
    "c5a_G": "#4E79A7",
    "bm213": "#0077BB",
}


def load_profile(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep=r"\s+", comment="#", names=["residue", "rmsf"])


def plot_profile(ax: plt.Axes, profiles: list[tuple[str, Path, str]], title: str, ylabel: str = "RMSF (Å)") -> None:
    for label, path, color_key in profiles:
        if not path.exists():
            continue
        df = load_profile(path)
        ax.plot(df["residue"], df["rmsf"], lw=1.7, color=PALETTE[color_key], label=label)
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("Residue number")
    ax.set_ylabel(ylabel)
    ax.legend(fontsize=7, frameon=False)
    style_axis(ax)


def draw_schematic(ax: plt.Axes) -> None:
    ax.set_xlim(0, 14)
    ax.set_ylim(0, 1)
    ax.axis("off")
    labels = ["H1", "ECL1", "H2", "ICL1", "H3", "ECL2", "H4", "ICL2", "H5", "ECL3", "H6", "ICL3", "H7", "H8"]
    colors = ["#F3F4F6", "#E5E7EB"] * 7
    for i, (label, color) in enumerate(zip(labels, colors)):
        ax.add_patch(Rectangle((i, 0.25), 1.0, 0.45, facecolor=color, edgecolor="black", lw=0.8))
        ax.text(i + 0.5, 0.475, label, ha="center", va="center", fontsize=8, fontweight="bold")
    ax.text(0.0, 0.95, "C5aR1 region map", ha="left", va="top", fontsize=11, fontweight="bold")
    ax.text(0.0, 0.08, "Schematic labels only; curves below use the available RMSF windows in this repo.", ha="left", va="bottom", fontsize=7)
    ax.text(-0.02, 0.47, "N-term", ha="right", va="center", fontsize=8)
    ax.text(14.02, 0.47, "C-term", ha="left", va="center", fontsize=8)


def main() -> None:
    ensure_dirs()
    apply_style()

    fig = plt.figure(figsize=(13.0, 8.8), constrained_layout=True)
    gs = fig.add_gridspec(3, 3, height_ratios=[0.62, 1.0, 1.0], hspace=0.42, wspace=0.28)

    ax = fig.add_subplot(gs[0, :])
    draw_schematic(ax)

    panel_specs = [
        ("a", "H3 / TM3 RMSF", [
            ("apo", Path("gpcr/apo/rmsf-TM3.dat"), "apo"),
            ("C5a", Path("gpcr/c5a/rmsf-TM3.dat"), "c5a"),
        ]),
        ("b", "TM6 / H6 RMSF", [
            ("apo", Path("gpcr/apo/rmsf-TM6.dat"), "apo"),
            ("C5a", Path("gpcr/c5a/rmsf-TM6.dat"), "c5a"),
        ]),
        ("c", "TM7 / H7 RMSF", [
            ("apo", Path("gpcr/apo/rmsf-TM7.dat"), "apo"),
            ("C5a", Path("gpcr/c5a/rmsf-TM7.dat"), "c5a"),
        ]),
        ("d", "ICL3 RMSF", [
            ("apo", Path("gpcr/apo/second/rmsf-icl3.dat"), "apo"),
            ("C5a_G", Path("gpcr/c5a_G/rmsf-icl3.dat"), "c5a_G"),
        ]),
        ("e", "H8 / hook RMSF", [
            ("C5a hook", Path("gpcr/c5a/rmsf-hook.dat"), "c5a"),
        ]),
        ("f", "TM4-TM5 core RMSF", [
            ("BM213 core", Path("gpcr/bm213/rmsf-protein.dat"), "bm213"),
        ]),
    ]

    for (label, title, profiles), grid_pos in zip(panel_specs, [gs[1, 0], gs[1, 1], gs[1, 2], gs[2, 0], gs[2, 1], gs[2, 2]]):
        ax = fig.add_subplot(grid_pos)
        ax.text(-0.08, 1.04, label, transform=ax.transAxes, ha="left", va="bottom", fontweight="bold", fontsize=12)
        plot_profile(ax, profiles, title)
        if label != "f":
            ax.set_xlabel("")

    fig.suptitle("C5aR1 segment-resolved RMSF profiles", fontsize=13, fontweight="bold")
    save_png(fig, BASIC_SUPP_DIR / "C5AR_rmsf_regions.png")


if __name__ == "__main__":
    main()
