#!/usr/bin/env python3
"""Run the C5AR beta-arrestin HADDOCK 9-cluster small screen."""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path


ROOT = Path("/mnt/e/work/modeling/c5ar/md")
DOCK_ROOT = ROOT / "article/analysis/docking/haddock"
FORMAL_DIR = ROOT / "article/analysis/sommd/pymol/formal"
TABLE_DIR = DOCK_ROOT / "tables"
CASE_DIR = DOCK_ROOT / "cases"

HADDOCK_HOME = Path("/mnt/e/app/haddock/haddock2.5-2026-03/haddock2.5-2026-03")
HADDOCK_PY = Path("/mnt/e/app/conda-envs/haddock2.5/bin/python")
RUN_HADDOCK = HADDOCK_HOME / "haddock/run_haddock.py"

RECEPTOR_ACTIVE = [108, 109, 207, 210, 271, 272]
LIGAND_ACTIVE = [67, 71, 73, 244, 249]
STRUCTURES_0 = 10
STRUCTURES_1 = 5
WATERREFINE = 5

HIS_ALIASES = {"HSD": "HIS", "HSE": "HIS", "HSP": "HIS", "HID": "HIS", "HIE": "HIS", "HIP": "HIS"}


@dataclass(frozen=True)
class ClusterCase:
    cluster: int
    system: str
    receptor: Path

    @property
    def name(self) -> str:
        return f"cluster_{self.cluster:02d}"

    @property
    def case_dir(self) -> Path:
        return CASE_DIR / self.name


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--jobs", type=int, default=3, help="Concurrent HADDOCK cases.")
    parser.add_argument("--only-cluster", action="append", default=[], help="Cluster id(s), e.g. 1 or 1,2,3.")
    parser.add_argument("--overwrite", action="store_true", help="Remove and rebuild selected case directories.")
    parser.add_argument("--setup-only", action="store_true", help="Prepare run1/run.cns but do not run HADDOCK.")
    parser.add_argument("--no-run", action="store_true", help="Alias for --setup-only.")
    parser.add_argument("--keep-hydrogens", action="store_true", help="Keep hydrogens in staged input PDBs.")
    return parser.parse_args()


def selected_clusters(raw: list[str]) -> set[int] | None:
    if not raw:
        return None
    out: set[int] = set()
    for item in raw:
        for part in item.split(","):
            part = part.strip()
            if part:
                out.add(int(part))
    return out


def discover_cases(only: set[int] | None = None) -> list[ClusterCase]:
    cases: list[ClusterCase] = []
    for pdb in sorted(FORMAL_DIR.glob("cluster_*.pdb")):
        match = re.match(r"cluster_(\d+)_([^_]+)_.*\.pdb$", pdb.name)
        if not match:
            continue
        cluster = int(match.group(1))
        if only is not None and cluster not in only:
            continue
        cases.append(ClusterCase(cluster=cluster, system=match.group(2), receptor=pdb))
    missing = sorted(set(range(1, 10)) - {c.cluster for c in cases}) if only is None else []
    if missing:
        raise FileNotFoundError(f"Missing formal cluster PDBs for clusters: {missing}")
    if not cases:
        raise FileNotFoundError("No selected cluster PDBs found.")
    return cases


def haddock_env() -> dict[str, str]:
    env = os.environ.copy()
    env["HADDOCK"] = str(HADDOCK_HOME)
    env["HADDOCKTOOLS"] = str(HADDOCK_HOME / "tools")
    env["PYTHONPATH"] = str(HADDOCK_HOME)
    env["PATH"] = f"{HADDOCK_PY.parent}:{env.get('PATH', '')}"
    return env


def is_hydrogen(line: str) -> bool:
    element = line[76:78].strip().upper() if len(line) >= 78 else ""
    atom = line[12:16].strip().upper()
    return element == "H" or atom.startswith("H") or (len(atom) > 1 and atom[0].isdigit() and atom[1] == "H")


def set_atom_name(line: str, atom_name: str) -> str:
    line = line.rstrip("\n").ljust(80)
    return f"{line[:12]}{atom_name:>4}{line[16:80]}\n"


def stage_pdb(src: Path, dst: Path, chain: str, keep_hydrogens: bool = False) -> dict[str, int]:
    n_in = 0
    n_out = 0
    residues: set[tuple[str, int]] = set()
    with src.open() as handle, dst.open("w") as out:
        out.write(f"REMARK staged from {src}\n")
        for raw in handle:
            if not raw.startswith(("ATOM  ", "HETATM")):
                continue
            n_in += 1
            if not keep_hydrogens and is_hydrogen(raw):
                continue
            line = raw.rstrip("\n").ljust(80)
            resname = line[17:20].strip().upper()
            resname = HIS_ALIASES.get(resname, resname)
            atom_name = line[12:16].strip().upper()
            if atom_name == "OT1":
                line = set_atom_name(line, "O").rstrip("\n").ljust(80)
            elif atom_name == "OT2":
                line = set_atom_name(line, "OXT").rstrip("\n").ljust(80)
            line = f"{line[:17]}{resname:>3}{line[20:]}"
            line = f"{line[:21]}{chain}{line[22:]}"
            line = f"{line[:72]}{chain:>4}{line[76:]}"
            out.write(line[:80] + "\n")
            resid_text = line[22:26].strip()
            if resid_text.lstrip("-").isdigit():
                residues.add((chain, int(resid_text)))
            n_out += 1
        out.write("TER\nEND\n")
    return {"atoms_in": n_in, "atoms_out": n_out, "residues": len(residues)}


def residue_set(pdb: Path) -> set[tuple[str, int]]:
    residues: set[tuple[str, int]] = set()
    with pdb.open() as handle:
        for line in handle:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            chain = line[21].strip()
            resid = line[22:26].strip()
            if chain and resid.lstrip("-").isdigit():
                residues.add((chain, int(resid)))
    return residues


def ensure_active_residues(receptor: Path, ligand: Path) -> None:
    receptor_res = residue_set(receptor)
    ligand_res = residue_set(ligand)
    missing_receptor = [r for r in RECEPTOR_ACTIVE if ("X", r) not in receptor_res]
    missing_ligand = [r for r in LIGAND_ACTIVE if ("A", r) not in ligand_res]
    if missing_receptor or missing_ligand:
        raise ValueError(
            f"Missing active residues in staged inputs. receptor={missing_receptor}; ligand={missing_ligand}"
        )


def write_air_tbl(path: Path) -> None:
    ligand_terms = [f"        ( resid {resid:<3d} and segid A)" for resid in LIGAND_ACTIVE]
    with path.open("w") as out:
        out.write("! C5AR beta-arrestin active-site ambiguous AIR restraints\n")
        out.write("! Receptor segid X active residues against beta-arrestin segid A active residues\n!\n")
        for rec_resid in RECEPTOR_ACTIVE:
            out.write(f"assign ( resid {rec_resid:<3d} and segid X)\n")
            out.write("       (\n")
            for idx, term in enumerate(ligand_terms):
                if idx:
                    out.write("     or\n")
                out.write(term + "\n")
            out.write("       )  2.0 2.0 0.0\n!\n")


def write_run_param(path: Path) -> None:
    text = "\n".join(
        [
            "AMBIG_TBL=./ambig.tbl",
            f"HADDOCK_DIR={HADDOCK_HOME}",
            "N_COMP=2",
            "PDB_FILE1=./receptor.pdb",
            "PDB_FILE2=./b-arrestin.pdb",
            "PROJECT_DIR=./",
            "PROT_SEGID_1=X",
            "PROT_SEGID_2=A",
            "RUN_NUMBER=1",
            "",
        ]
    )
    path.write_text(text)


def patch_run_cns(path: Path) -> None:
    replacements = {
        "structures_0": f"{{===>}} structures_0={STRUCTURES_0};",
        "structures_1": f"{{===>}} structures_1={STRUCTURES_1};",
        "waterrefine": f"{{===>}} waterrefine={WATERREFINE};",
        "firstwater": '{===>} firstwater="yes";',
    }
    seen = set()
    out_lines: list[str] = []
    for line in path.read_text().splitlines():
        stripped = line.strip()
        replaced = False
        for key, value in replacements.items():
            if re.match(rf"^\{{===>\}}\s*{re.escape(key)}=", stripped):
                out_lines.append(value)
                seen.add(key)
                replaced = True
                break
        if not replaced:
            out_lines.append(line)
    missing = sorted(set(replacements) - seen)
    if missing:
        raise ValueError(f"Could not patch run.cns parameters: {missing}")
    path.write_text("\n".join(out_lines) + "\n")


def run_command(cmd: list[str], cwd: Path, stdout_path: Path, stderr_path: Path) -> int:
    with stdout_path.open("w") as out, stderr_path.open("w") as err:
        proc = subprocess.run(cmd, cwd=cwd, env=haddock_env(), stdout=out, stderr=err, text=True)
    return proc.returncode


def water_models(case_dir: Path) -> list[Path]:
    water_dir = case_dir / "run1/structures/it1/water"
    file_list = water_dir / "file.list"
    models: list[Path] = []
    if file_list.exists():
        for raw in file_list.read_text(errors="ignore").splitlines():
            parts = raw.replace('"', "").split()
            if not parts:
                continue
            name = parts[0].replace("PREVIT:", "")
            candidate = water_dir / Path(name).name
            if candidate.exists() and candidate.suffix.lower() == ".pdb":
                models.append(candidate)
    if not models and water_dir.exists():
        models = sorted(p for p in water_dir.glob("*.pdb") if p.is_file())
    return models


def is_complete(case_dir: Path) -> bool:
    return (case_dir / "run1/structures/it1/water/file.list").exists() and len(water_models(case_dir)) > 0


def prepare_case(case: ClusterCase, overwrite: bool, keep_hydrogens: bool) -> dict[str, object]:
    case_dir = case.case_dir
    if overwrite and case_dir.exists():
        shutil.rmtree(case_dir)
    (case_dir / "logs").mkdir(parents=True, exist_ok=True)
    receptor_staged = case_dir / "receptor.pdb"
    ligand_staged = case_dir / "b-arrestin.pdb"
    ligand_src = FORMAL_DIR / "b-arrestin.pdb"
    if not ligand_src.exists():
        raise FileNotFoundError(ligand_src)
    receptor_stats = stage_pdb(case.receptor, receptor_staged, "X", keep_hydrogens=keep_hydrogens)
    ligand_stats = stage_pdb(ligand_src, ligand_staged, "A", keep_hydrogens=keep_hydrogens)
    ensure_active_residues(receptor_staged, ligand_staged)
    write_air_tbl(case_dir / "ambig.tbl")
    write_run_param(case_dir / "run.param")
    return {"receptor_stats": receptor_stats, "ligand_stats": ligand_stats}


def run_case(case: ClusterCase, overwrite: bool, setup_only: bool, keep_hydrogens: bool) -> dict[str, object]:
    start = time.time()
    row: dict[str, object] = {
        "cluster": case.cluster,
        "system": case.system,
        "receptor": str(case.receptor),
        "case_dir": str(case.case_dir),
        "status": "unknown",
        "n_water_models": 0,
        "elapsed_sec": 0.0,
        "setup_returncode": "",
        "run_returncode": "",
        "message": "",
    }
    try:
        if is_complete(case.case_dir) and not overwrite and not setup_only:
            row["status"] = "complete_existing"
            row["n_water_models"] = len(water_models(case.case_dir))
            return row

        stats = prepare_case(case, overwrite=overwrite, keep_hydrogens=keep_hydrogens)
        row.update(stats)
        run_cns = case.case_dir / "run1/run.cns"
        if not run_cns.exists():
            rc = run_command(
                [str(HADDOCK_PY), str(RUN_HADDOCK)],
                case.case_dir,
                case.case_dir / "logs/setup.stdout.log",
                case.case_dir / "logs/setup.stderr.log",
            )
            row["setup_returncode"] = rc
            if rc != 0 or not run_cns.exists():
                row["status"] = "setup_failed"
                row["message"] = "run.cns was not generated"
                return row
        else:
            row["setup_returncode"] = "existing"

        patch_run_cns(run_cns)
        if setup_only:
            row["status"] = "setup_complete"
            return row

        rc = run_command(
            [
                str(HADDOCK_PY),
                str(RUN_HADDOCK),
                "--batch_it0",
                str(STRUCTURES_0),
                "--batch_it1",
                str(STRUCTURES_1),
                "--batch_it1w",
                str(WATERREFINE),
            ],
            case.case_dir / "run1",
            case.case_dir / "logs/haddock.stdout.log",
            case.case_dir / "logs/haddock.stderr.log",
        )
        row["run_returncode"] = rc
        row["n_water_models"] = len(water_models(case.case_dir))
        if rc == 0 and is_complete(case.case_dir):
            row["status"] = "complete"
        elif is_complete(case.case_dir):
            row["status"] = "complete_with_nonzero_returncode"
            row["message"] = "water models exist but run_haddock returned nonzero"
        else:
            row["status"] = "run_failed"
            row["message"] = "water file.list/models not found"
    except Exception as exc:  # noqa: BLE001
        row["status"] = "error"
        row["message"] = f"{type(exc).__name__}: {exc}"
    finally:
        row["elapsed_sec"] = round(time.time() - start, 3)
    return row


def write_status(rows: list[dict[str, object]]) -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    fields = [
        "cluster",
        "system",
        "status",
        "n_water_models",
        "elapsed_sec",
        "setup_returncode",
        "run_returncode",
        "message",
        "case_dir",
        "receptor",
        "receptor_stats",
        "ligand_stats",
    ]
    with (TABLE_DIR / "haddock_batch_status.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in sorted(rows, key=lambda r: int(r["cluster"])):
            writer.writerow(row)


def write_summary_json(args: argparse.Namespace, cases: list[ClusterCase], rows: list[dict[str, object]]) -> None:
    summary = {
        "workflow": "C5AR beta-arrestin HADDOCK 9-cluster small screen",
        "constraint_mode": "active_site_ambiguous_air",
        "receptor_active": RECEPTOR_ACTIVE,
        "ligand_active": LIGAND_ACTIVE,
        "structures_0": STRUCTURES_0,
        "structures_1": STRUCTURES_1,
        "waterrefine": WATERREFINE,
        "firstwater": "yes",
        "jobs": args.jobs,
        "keep_hydrogens": bool(args.keep_hydrogens),
        "output_figure_format": "png_only",
        "haddock_home": str(HADDOCK_HOME),
        "cases": [case.name for case in cases],
        "status_counts": {status: sum(1 for r in rows if r["status"] == status) for status in sorted({r["status"] for r in rows})},
    }
    (DOCK_ROOT / "run-summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")


def main() -> int:
    args = parse_args()
    args.setup_only = args.setup_only or args.no_run
    if not RUN_HADDOCK.exists():
        raise FileNotFoundError(RUN_HADDOCK)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    CASE_DIR.mkdir(parents=True, exist_ok=True)
    cases = discover_cases(selected_clusters(args.only_cluster))
    rows: list[dict[str, object]] = []
    jobs = max(1, min(args.jobs, len(cases)))
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        futures = {
            pool.submit(
                run_case,
                case,
                args.overwrite,
                args.setup_only,
                args.keep_hydrogens,
            ): case
            for case in cases
        }
        for future in as_completed(futures):
            case = futures[future]
            row = future.result()
            rows.append(row)
            print(
                f"{case.name}: {row['status']} "
                f"models={row['n_water_models']} elapsed={row['elapsed_sec']}s",
                flush=True,
            )
            write_status(rows)
    write_status(rows)
    write_summary_json(args, cases, rows)
    ok_status = {"complete", "complete_existing", "complete_with_nonzero_returncode", "setup_complete"}
    return 0 if all(row["status"] in ok_status for row in rows) else 1


if __name__ == "__main__":
    sys.exit(main())
