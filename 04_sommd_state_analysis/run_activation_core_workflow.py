#!/usr/bin/env python3
from __future__ import annotations

import argparse
import subprocess
from pathlib import Path


ROOT = Path("/mnt/e/work/modeling/c5ar/md")
SOMMD_ROOT = ROOT / "article/analysis/sommd"
PYTHON = Path("/mnt/e/app/conda-envs/alloviz/bin/python")
R_PREFIX = Path("/mnt/e/app/conda-envs/R")

SCRIPT_DIR = SOMMD_ROOT / "scripts"
PREPARE = SCRIPT_DIR / "prepare_activation_core_features.py"
RUN_R = SCRIPT_DIR / "run_activation_core_som.R"
SUMMARY = SCRIPT_DIR / "summarize_activation_core_results.py"
OFFICIAL_TRANSITION = SCRIPT_DIR / "plot_official_transition_network.R"
FORMAL_POPULATION_FIGURES = SCRIPT_DIR / "plot_formal_som_population_figures.R"

DEFAULT_RLEN = {
    "smoke": 200,
    "pilot": 500,
    "formal": 1000,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run C5aR1 activation-core SOMMD workflow.")
    parser.add_argument("--profile", choices=("smoke", "pilot", "formal", "all"), default="smoke")
    parser.add_argument("--grid", default="8x8")
    parser.add_argument("--clusters", type=int, default=8)
    parser.add_argument("--seed", type=int, default=20260603)
    parser.add_argument("--toroidal", choices=("true", "false"), default="true")
    parser.add_argument("--scale", choices=("true", "false"), default="true")
    parser.add_argument("--rlen", type=int, default=None, help="Override default rlen for selected profile(s).")
    parser.add_argument("--force", action="store_true", help="Rebuild feature matrices.")
    return parser.parse_args()


def run(cmd: list[str]) -> None:
    print("[run]", " ".join(cmd), flush=True)
    subprocess.run(cmd, cwd=ROOT, check=True)


def run_profile(profile: str, args: argparse.Namespace) -> None:
    input_dir = SOMMD_ROOT / "inputs" / profile / "activation_core"
    run_dir = SOMMD_ROOT / "runs" / profile / "activation_core"
    rlen = args.rlen if args.rlen is not None else DEFAULT_RLEN[profile]

    prepare_cmd = [
        str(PYTHON),
        str(PREPARE),
        "--profile",
        profile,
    ]
    if args.force:
        prepare_cmd.append("--force")
    run(prepare_cmd)

    run(
        [
            "conda",
            "run",
            "-p",
            str(R_PREFIX),
            "Rscript",
            str(RUN_R),
            "--feature-matrix",
            str(input_dir / "feature_matrix.csv"),
            "--frame-metadata",
            str(input_dir / "frame_metadata.csv"),
            "--feature-metadata",
            str(input_dir / "feature_metadata.csv"),
            "--out",
            str(run_dir),
            "--grid",
            args.grid,
            "--clusters",
            str(args.clusters),
            "--rlen",
            str(rlen),
            "--seed",
            str(args.seed),
            "--toroidal",
            args.toroidal,
            "--scale",
            args.scale,
        ]
    )

    run(
        [
            str(PYTHON),
            str(SUMMARY),
            "--profile",
            profile,
            "--run-dir",
            str(run_dir),
            "--input-dir",
            str(input_dir),
        ]
    )

    if profile == "formal":
        run(
            [
                "conda",
                "run",
                "-p",
                str(R_PREFIX),
                "Rscript",
                str(OFFICIAL_TRANSITION),
                "--profile",
                profile,
                "--run_dir",
                str(run_dir),
                "--out_figures",
                str(SOMMD_ROOT / "figures"),
                "--out_tables",
                str(SOMMD_ROOT / "tables"),
                "--seed",
                str(args.seed),
            ]
        )
        run(
            [
                "conda",
                "run",
                "-p",
                str(R_PREFIX),
                "Rscript",
                str(FORMAL_POPULATION_FIGURES),
            ]
        )


def main() -> int:
    args = parse_args()
    profiles = ("smoke", "pilot", "formal") if args.profile == "all" else (args.profile,)
    for profile in profiles:
        run_profile(profile, args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
