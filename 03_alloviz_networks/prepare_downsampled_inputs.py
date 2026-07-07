#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import MDAnalysis as mda


ROOT = Path("/mnt/e/work/modeling/c5ar/md")
OUTROOT = ROOT / "article/analysis/alloviz"


@dataclass(frozen=True)
class SystemSpec:
    key: str
    pdb: str
    trajs: tuple[str, ...]


SYSTEMS = [
    SystemSpec(
        key="apo",
        pdb="gpcr/apo/md-initial.pdb",
        trajs=(
            "gpcr/apo/2000ns.xtc",
            "gpcr/apo/second/2000ns.xtc",
            "gpcr/apo/third/2000ns.xtc",
        ),
    ),
    SystemSpec(
        key="BM213",
        pdb="gpcr/bm213/md-initial.pdb",
        trajs=(
            "gpcr/bm213/2wframe.xtc",
            "gpcr/bm213/second/2wframe.xtc",
            "gpcr/bm213/third/2wframe.xtc",
        ),
    ),
    SystemSpec(
        key="C5apep",
        pdb="gpcr/c5apep/md-initial.pdb",
        trajs=(
            "gpcr/c5apep/2wframe.xtc",
            "gpcr/c5apep/second/2wframe2.xtc",
            "gpcr/c5apep/third/2wframe3.xtc",
        ),
    ),
]


def write_downsampled(profile: str, stride: int, force: bool) -> None:
    outdir = OUTROOT / "inputs/downsampled" / profile
    outdir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []

    for spec in SYSTEMS:
        pdb = ROOT / spec.pdb
        for rep_idx, traj_rel in enumerate(spec.trajs, start=1):
            traj = ROOT / traj_rel
            out = outdir / f"{spec.key}_rep{rep_idx}_stride{stride}.xtc"
            u = mda.Universe(str(pdb), str(traj))
            source_frames = len(u.trajectory)
            expected_frames = len(range(0, source_frames, stride))
            if out.exists() and not force:
                produced_frames = len(mda.Universe(str(pdb), str(out)).trajectory)
            else:
                with mda.Writer(str(out), n_atoms=u.atoms.n_atoms) as writer:
                    for ts in u.trajectory[::stride]:
                        writer.write(u.atoms)
                produced_frames = len(mda.Universe(str(pdb), str(out)).trajectory)

            rows.append(
                {
                    "profile": profile,
                    "system": spec.key,
                    "replicate": f"rep{rep_idx}",
                    "pdb": spec.pdb,
                    "source_trajectory": traj_rel,
                    "downsampled_trajectory": out.relative_to(ROOT).as_posix(),
                    "stride": stride,
                    "source_frames": source_frames,
                    "expected_frames": expected_frames,
                    "produced_frames": produced_frames,
                }
            )

    manifest = outdir / "downsample_manifest.csv"
    with manifest.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(manifest)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--profile", choices=["smoke", "pilot", "formal"], default="pilot")
    parser.add_argument("--stride", type=int, default=None)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    default_stride = {"smoke": 500, "pilot": 100, "formal": 10}[args.profile]
    write_downsampled(args.profile, args.stride or default_stride, args.force)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
