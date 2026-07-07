#!/usr/bin/env python3
"""Validate the C5aR 2 us trajectory manifest.

The script checks path existence, atom-count compatibility, frame counts,
trajectory spacing, approximate sampled time, and basic receptor/ligand
selection sizes. It does not modify source trajectories.
"""

from __future__ import annotations

import argparse
import csv
import sys
import warnings
from pathlib import Path

import MDAnalysis as mda


DEFAULT_MANIFEST = Path("article/analysis/c5ar_2us_inputs/trajectory_manifest.tsv")
DEFAULT_OUTPUT = Path("article/analysis/c5ar_2us_inputs/input_validation.tsv")


def as_float(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def validate_row(row: dict[str, str], root: Path) -> dict[str, str]:
    structure = root / row["structure"]
    trajectory = root / row["trajectory"]
    topology = root / row["topology"]

    result = {
        "system": row["system"],
        "replicate": row["replicate"],
        "structure": row["structure"],
        "trajectory": row["trajectory"],
        "structure_exists": str(structure.is_file()),
        "trajectory_exists": str(trajectory.is_file()),
        "topology_exists": str(topology.is_file()),
        "structure_atoms": "NA",
        "trajectory_atoms": "NA",
        "atom_count_match": "NA",
        "frames": "NA",
        "dt_ps": "NA",
        "sampled_ns_frames_x_dt": "NA",
        "last_time_ns": "NA",
        "receptor_atoms": "NA",
        "ligand_atoms": "NA",
        "status": "unchecked",
        "message": "",
    }

    if not structure.is_file() or not trajectory.is_file():
        result["status"] = "missing_input"
        result["message"] = "Structure or trajectory path is missing."
        return result

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            universe = mda.Universe(str(structure), str(trajectory))
            structure_atoms = universe.atoms.n_atoms
            trajectory_atoms = universe.trajectory.n_atoms
            frames = len(universe.trajectory)
            dt_ps = float(universe.trajectory.dt)
            universe.trajectory[-1]
            last_time_ns = float(universe.trajectory.time) / 1000.0

        receptor_atoms = universe.select_atoms(row["receptor_selection"]).n_atoms
        ligand_atoms = universe.select_atoms(row["ligand_selection"]).n_atoms

        result.update(
            {
                "structure_atoms": str(structure_atoms),
                "trajectory_atoms": str(trajectory_atoms),
                "atom_count_match": str(structure_atoms == trajectory_atoms),
                "frames": str(frames),
                "dt_ps": f"{dt_ps:.6g}",
                "sampled_ns_frames_x_dt": f"{frames * dt_ps / 1000.0:.6g}",
                "last_time_ns": f"{last_time_ns:.6g}",
                "receptor_atoms": str(receptor_atoms),
                "ligand_atoms": str(ligand_atoms),
            }
        )

        expected_frames = int(row["source_frames"])
        expected_dt = as_float(row["source_dt_ps"])
        intended_ns = as_float(row["intended_ns"])
        sampled_ns = frames * dt_ps / 1000.0

        checks = [
            structure_atoms == trajectory_atoms,
            frames == expected_frames,
            abs(dt_ps - expected_dt) < 1e-3,
            abs(sampled_ns - intended_ns) <= 0.2,
            receptor_atoms > 0,
            ligand_atoms > 0,
        ]
        result["status"] = "ok" if all(checks) else "check"
        if result["status"] == "check":
            failed = []
            if structure_atoms != trajectory_atoms:
                failed.append("atom_count")
            if frames != expected_frames:
                failed.append("frames")
            if abs(dt_ps - expected_dt) >= 1e-3:
                failed.append("dt")
            if abs(sampled_ns - intended_ns) > 0.2:
                failed.append("sampled_ns")
            if receptor_atoms <= 0:
                failed.append("receptor_selection")
            if ligand_atoms <= 0:
                failed.append("ligand_selection")
            result["message"] = ",".join(failed)
    except Exception as exc:  # noqa: BLE001 - report row-level validation failures.
        result["status"] = "error"
        result["message"] = f"{type(exc).__name__}: {exc}"

    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    root = Path.cwd()
    if not args.manifest.is_file():
        print(f"Manifest not found: {args.manifest}", file=sys.stderr)
        return 2

    with args.manifest.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    results = [validate_row(row, root) for row in rows]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(results[0].keys()),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(results)

    status_counts: dict[str, int] = {}
    for result in results:
        status_counts[result["status"]] = status_counts.get(result["status"], 0) + 1

    print(f"Wrote {args.output}")
    for status, count in sorted(status_counts.items()):
        print(f"{status}: {count}")

    return 0 if all(result["status"] == "ok" for result in results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
