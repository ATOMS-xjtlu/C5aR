#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import MDAnalysis as mda


ROOT = Path("/mnt/e/work/modeling/c5ar/md")
OUTROOT = ROOT / "article/analysis/alloviz"
PYTHON = Path("/mnt/e/app/conda-envs/alloviz/bin/python")
ALLOVIZ_RUN = OUTROOT / "scripts/alloviz_run_serial.py"


@dataclass(frozen=True)
class SystemSpec:
    key: str
    label: str
    pdb: str
    trajs: tuple[str, ...]
    protein_sel: str
    role: str


@dataclass(frozen=True)
class RunSpec:
    system: SystemSpec
    replicate: str
    traj: str


SYSTEMS = [
    SystemSpec(
        key="apo",
        label="apo",
        pdb="gpcr/apo/md-initial.pdb",
        trajs=(
            "gpcr/apo/2000ns.xtc",
            "gpcr/apo/second/2000ns.xtc",
            "gpcr/apo/third/2000ns.xtc",
        ),
        protein_sel="protein",
        role="basal reference",
    ),
    SystemSpec(
        key="BM213",
        label="BM213",
        pdb="gpcr/bm213/md-initial.pdb",
        trajs=(
            "gpcr/bm213/2wframe.xtc",
            "gpcr/bm213/second/2wframe.xtc",
            "gpcr/bm213/third/2wframe.xtc",
        ),
        protein_sel="protein and not resid 282-289",
        role="G protein-biased / stronger agonist state",
    ),
    SystemSpec(
        key="C5apep",
        label="C5apep",
        pdb="gpcr/c5apep/md-initial.pdb",
        trajs=(
            "gpcr/c5apep/2wframe.xtc",
            "gpcr/c5apep/second/2wframe2.xtc",
            "gpcr/c5apep/third/2wframe3.xtc",
        ),
        protein_sel="protein",
        role="partial efficacy state",
    ),
]

PILOT_PKGS = ("GetContacts", "correlationplus_CA_Pear", "pytraj_CA")
FORMAL_REQUIRED = ("GetContacts",)
FORMAL_OPTIONAL = ("correlationplus_CA_Pear", "pytraj_CA")


def iter_runs() -> list[RunSpec]:
    runs: list[RunSpec] = []
    for spec in SYSTEMS:
        for idx, traj in enumerate(spec.trajs, start=1):
            runs.append(RunSpec(system=spec, replicate=f"rep{idx}", traj=traj))
    return runs


def downsampled_traj(run: RunSpec, phase: str) -> str:
    profile = "pilot" if phase == "pilot" else "formal"
    base = OUTROOT / "inputs/downsampled" / profile
    matches = sorted(base.glob(f"{run.system.key}_{run.replicate}_stride*.xtc"))
    if matches:
        return matches[0].relative_to(ROOT).as_posix()
    return run.traj


def input_tag(run: RunSpec, phase: str) -> str:
    base = "downsampled" if downsampled_traj(run, phase) != run.traj else "raw"
    return f"serial_c1_{base}"


def run_outdir(run: RunSpec, phase: str, pkg: str) -> Path:
    return OUTROOT / "runs" / phase / run.system.key / run.replicate / f"{pkg}_{input_tag(run, phase)}"


def ensure_dirs() -> None:
    for rel in ("inputs", "runs", "tables", "figures", "logs"):
        (OUTROOT / rel).mkdir(parents=True, exist_ok=True)


def validate_systems() -> None:
    manifest_rows: list[dict[str, object]] = []
    receptor_rows: list[dict[str, object]] = []

    for spec in SYSTEMS:
        pdb = ROOT / spec.pdb
        if not pdb.exists():
            raise FileNotFoundError(pdb)

        u_top = mda.Universe(str(pdb))
        receptor = u_top.select_atoms(spec.protein_sel)
        if receptor.n_atoms == 0:
            raise ValueError(f"{spec.key}: receptor selection returned zero atoms: {spec.protein_sel}")

        for residue_idx, residue in enumerate(receptor.residues, start=1):
            receptor_rows.append(
                {
                    "system": spec.key,
                    "raw_resid": residue.resid,
                    "resname": residue.resname,
                    "residue_index_1based": residue_idx,
                }
            )

        for rep_idx, traj_rel in enumerate(spec.trajs, start=1):
            traj = ROOT / traj_rel
            if not traj.exists():
                raise FileNotFoundError(traj)

            u = mda.Universe(str(pdb), str(traj))
            first_time = float(u.trajectory[0].time)
            last_time = float(u.trajectory[-1].time)
            selected = u.select_atoms(spec.protein_sel)
            if u.atoms.n_atoms != u_top.atoms.n_atoms:
                raise ValueError(
                    f"{spec.key} rep{rep_idx}: topology atoms {u_top.atoms.n_atoms} "
                    f"!= trajectory atoms {u.atoms.n_atoms}"
                )
            if selected.n_atoms != receptor.n_atoms:
                raise ValueError(
                    f"{spec.key} rep{rep_idx}: selected atoms {selected.n_atoms} "
                    f"!= receptor atoms {receptor.n_atoms}"
                )

            manifest_rows.append(
                {
                    "system": spec.key,
                    "label": spec.label,
                    "role": spec.role,
                    "replicate": f"rep{rep_idx}",
                    "pdb": spec.pdb,
                    "trajectory": traj_rel,
                    "protein_sel": spec.protein_sel,
                    "topology_atoms": u_top.atoms.n_atoms,
                    "trajectory_atoms": u.atoms.n_atoms,
                    "receptor_atoms": receptor.n_atoms,
                    "receptor_residues": len(receptor.residues),
                    "frames": len(u.trajectory),
                    "first_time_ps": first_time,
                    "last_time_ps": last_time,
                }
            )

    manifest_csv = OUTROOT / "inputs/input_manifest.csv"
    with manifest_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(manifest_rows[0].keys()))
        writer.writeheader()
        writer.writerows(manifest_rows)

    receptor_csv = OUTROOT / "inputs/receptor_selection_residues.csv"
    with receptor_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(receptor_rows[0].keys()))
        writer.writeheader()
        writer.writerows(receptor_rows)

    manifest_json = OUTROOT / "inputs/run_config.json"
    manifest_json.write_text(
        json.dumps(
            {
                "created_at": datetime.now().isoformat(timespec="seconds"),
                "alloviz_python": str(PYTHON),
                "alloviz_helper": str(ALLOVIZ_RUN),
                "pilot_pkgs": PILOT_PKGS,
                "formal_required": FORMAL_REQUIRED,
                "formal_optional": FORMAL_OPTIONAL,
                "systems": [
                    {
                        "key": spec.key,
                        "label": spec.label,
                        "role": spec.role,
                        "pdb": spec.pdb,
                        "trajs": list(spec.trajs),
                        "protein_sel": spec.protein_sel,
                    }
                    for spec in SYSTEMS
                ],
            },
            indent=2,
            ensure_ascii=False,
        )
        + "\n",
        encoding="utf-8",
    )


def export_exists(outdir: Path, pkg: str) -> bool:
    return (
        (outdir / "run_summary.json").exists()
        and (outdir / "exports" / f"{pkg}__Spatially_distant__nodes.csv").exists()
        and (outdir / "exports" / f"{pkg}__Spatially_distant__edges.csv").exists()
    )


def run_one(run: RunSpec, phase: str, pkg: str, stride: int | None, force: bool) -> dict[str, object]:
    spec = run.system
    traj_for_run = downsampled_traj(run, phase)
    outdir = run_outdir(run, phase, pkg)
    outdir.mkdir(parents=True, exist_ok=True)
    log_prefix = OUTROOT / "logs" / f"{phase}_{spec.key}_{run.replicate}_{pkg}"
    status_path = outdir / "status.json"

    if not force and export_exists(outdir, pkg):
        status = {
            "system": spec.key,
            "replicate": run.replicate,
            "phase": phase,
            "pkg": pkg,
            "status": "skipped_existing",
            "outdir": str(outdir),
        }
        status_path.write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
        return status

    cmd = [
        str(PYTHON),
        str(ALLOVIZ_RUN),
        "--pdb",
        spec.pdb,
        "--trajs",
        traj_for_run,
        "--out",
        str(outdir),
        "--name",
        f"{spec.key}_{run.replicate}_{phase}_{pkg}",
        "--pkgs",
        pkg,
        "--filterings",
        "Spatially_distant",
        "--elements",
        "nodes,edges",
        "--metrics",
        "btw",
        "--cores",
        "1",
        "--taskcpus",
        "1",
        "--protein-sel",
        spec.protein_sel,
    ]
    if stride is not None and traj_for_run == run.traj:
        cmd.extend(["--stride", str(stride)])

    with (log_prefix.with_suffix(".out")).open("w", encoding="utf-8") as stdout, (
        log_prefix.with_suffix(".err")
    ).open("w", encoding="utf-8") as stderr:
        result = subprocess.run(cmd, cwd=str(ROOT), stdout=stdout, stderr=stderr, text=True)

    status = {
        "system": spec.key,
        "replicate": run.replicate,
        "phase": phase,
        "pkg": pkg,
        "stride": stride,
        "source_trajectory": run.traj,
        "trajectory": traj_for_run,
        "status": "success" if result.returncode == 0 and export_exists(outdir, pkg) else "failed",
        "returncode": result.returncode,
        "outdir": str(outdir),
        "stdout": str(log_prefix.with_suffix(".out")),
        "stderr": str(log_prefix.with_suffix(".err")),
        "command": cmd,
    }
    status_path.write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
    return status


def load_status(phase: str, run: RunSpec, pkg: str) -> dict[str, object] | None:
    path = run_outdir(run, phase, pkg) / "status.json"
    if not path.exists():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def is_success(status: dict[str, object] | None) -> bool:
    return bool(status and status.get("status") in {"success", "skipped_existing"})


def write_method_failures(statuses: list[dict[str, object]]) -> None:
    failures = [s for s in statuses if s.get("status") == "failed"]
    lines = ["# AlloViz method failures", ""]
    if not failures:
        lines.append("No AlloViz method failures were recorded.")
    else:
        lines.append("| phase | system | replicate | package | returncode | stderr |")
        lines.append("|---|---|---|---:|---:|---|")
        for item in failures:
            lines.append(
                f"| {item['phase']} | {item['system']} | {item.get('replicate', 'NA')} | {item['pkg']} | "
                f"{item.get('returncode', 'NA')} | `{item.get('stderr', 'NA')}` |"
            )
    (OUTROOT / "logs/method_failures.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase", choices=["pilot", "formal", "all"], default="all")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    ensure_dirs()
    validate_systems()

    statuses: list[dict[str, object]] = []
    if args.phase in {"pilot", "all"}:
        for run in iter_runs():
            for pkg in PILOT_PKGS:
                print(f"[pilot] {run.system.key} {run.replicate} {pkg}", flush=True)
                statuses.append(run_one(run, "pilot", pkg, stride=10, force=args.force))

    if args.phase in {"formal", "all"}:
        for run in iter_runs():
            for pkg in FORMAL_REQUIRED:
                print(f"[formal] {run.system.key} {run.replicate} {pkg}", flush=True)
                statuses.append(run_one(run, "formal", pkg, stride=None, force=args.force))
            for pkg in FORMAL_OPTIONAL:
                pilot_status = load_status("pilot", run, pkg)
                if args.phase == "formal" and pilot_status is None:
                    print(f"[formal-skip] {run.system.key} {run.replicate} {pkg}: no pilot status", flush=True)
                    continue
                if is_success(pilot_status):
                    print(f"[formal] {run.system.key} {run.replicate} {pkg}", flush=True)
                    statuses.append(run_one(run, "formal", pkg, stride=None, force=args.force))
                else:
                    status = {
                        "system": run.system.key,
                        "replicate": run.replicate,
                        "phase": "formal",
                        "pkg": pkg,
                        "status": "skipped_failed_pilot",
                        "pilot_status": pilot_status,
                    }
                    run_outdir(run, "formal", pkg).mkdir(parents=True, exist_ok=True)
                    (run_outdir(run, "formal", pkg) / "status.json").write_text(
                        json.dumps(status, indent=2) + "\n", encoding="utf-8"
                    )
                    statuses.append(status)

    expected_status_files: list[Path] = []
    phases = []
    if args.phase in {"pilot", "all"}:
        phases.append("pilot")
    if args.phase in {"formal", "all"}:
        phases.append("formal")
    for phase in phases:
        pkgs = PILOT_PKGS if phase == "pilot" else FORMAL_REQUIRED + FORMAL_OPTIONAL
        for run in iter_runs():
            for pkg in pkgs:
                expected_status_files.append(run_outdir(run, phase, pkg) / "status.json")

    all_statuses = []
    for status_file in expected_status_files:
        if status_file.exists():
            all_statuses.append(json.loads(status_file.read_text(encoding="utf-8")))
    write_method_failures(all_statuses)

    summary_path = OUTROOT / "logs/run_status_summary.json"
    summary_path.write_text(json.dumps(all_statuses, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {summary_path}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
