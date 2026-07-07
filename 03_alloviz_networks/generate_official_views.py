#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import traceback
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import AlloViz
from alloviz_run_serial import install_serial_analysis_patch


ROOT = Path("/mnt/e/work/modeling/c5ar/md")
OUTROOT = ROOT / "article/analysis/alloviz"
VIEWROOT = OUTROOT / "official_views"


@dataclass(frozen=True)
class ViewSystem:
    key: str
    pdb: str
    trajs: tuple[str, str, str]
    protein_sel: str


SYSTEMS = {
    "apo": ViewSystem(
        key="apo",
        pdb="gpcr/apo/md-initial.pdb",
        trajs=(
            "gpcr/apo/2000ns.xtc",
            "gpcr/apo/second/2000ns.xtc",
            "gpcr/apo/third/2000ns.xtc",
        ),
        protein_sel="protein",
    ),
    "BM213": ViewSystem(
        key="BM213",
        pdb="gpcr/bm213/md-initial.pdb",
        trajs=(
            "gpcr/bm213/2wframe.xtc",
            "gpcr/bm213/second/2wframe.xtc",
            "gpcr/bm213/third/2wframe.xtc",
        ),
        protein_sel="protein and not resid 282-289",
    ),
    "C5apep": ViewSystem(
        key="C5apep",
        pdb="gpcr/c5apep/md-initial.pdb",
        trajs=(
            "gpcr/c5apep/2wframe.xtc",
            "gpcr/c5apep/second/2wframe2.xtc",
            "gpcr/c5apep/third/2wframe3.xtc",
        ),
        protein_sel="protein",
    ),
}

SYSTEM_VIEW_PKG = "PyInteraph2_Contacts"
SYSTEM_VIEW_METRIC = "btw"
SYSTEM_VIEW_FILTERING = "Spatially_distant"
SYSTEM_VIEW_ELEMENT = "nodes"

DELTA_PKG = "correlationplus_CA_Pear"
DELTA_FILTERING = "No_Sequence_Neighbors"
DELTA_ELEMENT = "edges"
DELTA_METRIC = "btw"


def rel(path: str) -> str:
    return str((ROOT / path).resolve())


def preferred_traj(system: ViewSystem, replicate: int) -> str:
    rep = f"rep{replicate}"
    for profile in ("formal", "pilot", "smoke"):
        base = OUTROOT / "inputs/downsampled" / profile
        matches = sorted(base.glob(f"{system.key}_{rep}_stride*.xtc"))
        if matches:
            return matches[0].relative_to(ROOT).as_posix()
    return system.trajs[replicate - 1]


def create_protein(system: ViewSystem, replicate: int, outdir: Path, pkgs: list[str], filterings: list[str]) -> AlloViz.Protein:
    traj = preferred_traj(system, replicate)
    obj = AlloViz.Protein(
        pdb=rel(system.pdb),
        trajs=rel(traj),
        name=f"{system.key}_rep{replicate}_official_view",
        path=str(outdir),
        protein_sel=system.protein_sel,
    )
    obj.calculate(pkgs=pkgs, cores=1, taskcpus=1, stride=10)
    for filtering in filterings:
        obj.filter(
            pkgs=pkgs,
            filterings=filtering,
            GetContacts_threshold=0.0,
            Sequence_Neighbor_distance=5,
            Interresidue_distance=10.0,
        )
    obj.analyze(pkgs=pkgs, filterings=filterings, elements=["nodes", "edges"], metrics=DELTA_METRIC, cores=1)
    return obj


def safe_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def run_system_view(system_key: str, replicate: int, force: bool) -> dict[str, object]:
    system = SYSTEMS[system_key]
    outdir = VIEWROOT / "systems" / system_key / f"rep{replicate}"
    status_path = outdir / "view_status.json"
    if status_path.exists() and not force:
        return json.loads(status_path.read_text(encoding="utf-8"))

    status: dict[str, object] = {
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "type": "system.view",
        "system": system_key,
        "replicate": f"rep{replicate}",
        "pkg": SYSTEM_VIEW_PKG,
        "metric": SYSTEM_VIEW_METRIC,
        "filtering": SYSTEM_VIEW_FILTERING,
        "element": SYSTEM_VIEW_ELEMENT,
        "status": "started",
        "outdir": str(outdir),
    }
    try:
        protein = create_protein(
            system,
            replicate,
            outdir,
            pkgs=[SYSTEM_VIEW_PKG],
            filterings=[SYSTEM_VIEW_FILTERING],
        )
        view_result = protein.view(SYSTEM_VIEW_PKG, SYSTEM_VIEW_METRIC, SYSTEM_VIEW_FILTERING, SYSTEM_VIEW_ELEMENT)
        status.update({"status": "success", "view_result_repr": repr(view_result)})
    except Exception as exc:
        status.update(
            {
                "status": "failed",
                "error": f"{type(exc).__name__}: {exc}",
                "traceback": traceback.format_exc(),
            }
        )
    safe_write(status_path, json.dumps(status, indent=2, ensure_ascii=False) + "\n")
    return status


def run_delta_view(left_key: str, right_key: str, replicate: int, force: bool) -> dict[str, object]:
    outdir = VIEWROOT / "delta" / f"{left_key}_vs_{right_key}" / f"rep{replicate}"
    status_path = outdir / "delta_view_status.json"
    if status_path.exists() and not force:
        return json.loads(status_path.read_text(encoding="utf-8"))

    status: dict[str, object] = {
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "type": "AlloViz.Delta.view",
        "comparison": f"{left_key}_vs_{right_key}",
        "left": left_key,
        "right": right_key,
        "replicate": f"rep{replicate}",
        "pkg": DELTA_PKG,
        "filtering": DELTA_FILTERING,
        "element": DELTA_ELEMENT,
        "metric": DELTA_METRIC,
        "status": "started",
        "outdir": str(outdir),
    }
    try:
        left = create_protein(
            SYSTEMS[left_key],
            replicate,
            outdir / left_key,
            pkgs=[DELTA_PKG],
            filterings=[DELTA_FILTERING],
        )
        right = create_protein(
            SYSTEMS[right_key],
            replicate,
            outdir / right_key,
            pkgs=[DELTA_PKG],
            filterings=[DELTA_FILTERING],
        )
        delta = AlloViz.Delta(left, right)
        if hasattr(delta, "_aln"):
            safe_write(outdir / "delta_alignment.txt", str(delta._aln) + "\n")
        view_result = delta.view(pkg=DELTA_PKG, metric=DELTA_METRIC, filtering=DELTA_FILTERING, element=DELTA_ELEMENT)
        status.update({"status": "success", "view_result_repr": repr(view_result)})
    except Exception as exc:
        status.update(
            {
                "status": "failed",
                "error": f"{type(exc).__name__}: {exc}",
                "traceback": traceback.format_exc(),
            }
        )
    safe_write(status_path, json.dumps(status, indent=2, ensure_ascii=False) + "\n")
    return status


def write_manifest(statuses: list[dict[str, object]]) -> None:
    manifest = {
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "note": "Official AlloViz tutorial-style views. These are visualization outputs; statistical comparisons remain in tables/.",
        "system_view": {
            "call": 'system.view("PyInteraph2_Contacts", "btw", "Spatially_distant", "nodes")',
            "pkg": SYSTEM_VIEW_PKG,
            "metric": SYSTEM_VIEW_METRIC,
            "filtering": SYSTEM_VIEW_FILTERING,
            "element": SYSTEM_VIEW_ELEMENT,
        },
        "delta_view": {
            "call": 'delta.view(pkg="correlationplus_CA_Pear", metric="btw", filtering="No_Sequence_Neighbors", element="edges")',
            "pkg": DELTA_PKG,
            "filtering": DELTA_FILTERING,
            "element": DELTA_ELEMENT,
            "metric": DELTA_METRIC,
        },
        "statuses": statuses,
    }
    safe_write(VIEWROOT / "official_view_manifest.json", json.dumps(manifest, indent=2, ensure_ascii=False) + "\n")


def main() -> int:
    install_serial_analysis_patch()

    parser = argparse.ArgumentParser()
    parser.add_argument("--replicate", type=int, default=1, choices=[1, 2, 3])
    parser.add_argument("--all-repeats", action="store_true")
    parser.add_argument("--systems", action="store_true", help="Generate system.view outputs.")
    parser.add_argument("--deltas", action="store_true", help="Generate Delta.view outputs.")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    do_systems = args.systems or not args.deltas
    do_deltas = args.deltas or not args.systems
    replicates = [1, 2, 3] if args.all_repeats else [args.replicate]
    VIEWROOT.mkdir(parents=True, exist_ok=True)

    statuses: list[dict[str, object]] = []
    for rep in replicates:
        if do_systems:
            for system_key in ("apo", "BM213", "C5apep"):
                print(f"[system.view] {system_key} rep{rep}", flush=True)
                statuses.append(run_system_view(system_key, rep, force=args.force))
        if do_deltas:
            for left, right in (("BM213", "apo"), ("C5apep", "apo"), ("BM213", "C5apep")):
                print(f"[Delta.view] {left} vs {right} rep{rep}", flush=True)
                statuses.append(run_delta_view(left, right, rep, force=args.force))
    write_manifest(statuses)
    print(f"Wrote {VIEWROOT / 'official_view_manifest.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
