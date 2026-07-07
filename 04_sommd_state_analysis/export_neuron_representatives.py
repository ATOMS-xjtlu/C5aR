#!/usr/bin/env python3
from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

import MDAnalysis as mda
import pandas as pd


ROOT = Path("/mnt/e/work/modeling/c5ar/md")
SOMMD_ROOT = ROOT / "article/analysis/sommd"
R_PREFIX = Path("/mnt/e/app/conda-envs/R")
SYSTEM_COLORS = {"apo": "gray70", "BM213": "tv_blue", "C5apep": "orange"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Export one representative PDB for each occupied SOM neuron.")
    parser.add_argument("--profile", default="formal", choices=("smoke", "pilot", "formal"))
    parser.add_argument("--root", default=str(ROOT))
    parser.add_argument("--r-prefix", default=str(R_PREFIX))
    parser.add_argument("--skip-select", action="store_true", help="Reuse the existing neuron representative CSV.")
    return parser.parse_args()


def run_selection(root: Path, profile: str, r_prefix: Path, out_csv: Path) -> None:
    script = root / "article/analysis/sommd/scripts/select_neuron_representatives.R"
    cmd = [
        "conda",
        "run",
        "-p",
        str(r_prefix),
        "Rscript",
        str(script),
        "--profile",
        profile,
        "--root",
        str(root),
        "--out-csv",
        str(out_csv),
    ]
    subprocess.run(cmd, check=True)


def relative(path: Path, root: Path) -> str:
    return path.relative_to(root).as_posix()


def export_pdbs(root: Path, profile: str) -> tuple[Path, Path, int]:
    sommd_root = root / "article/analysis/sommd"
    table_path = sommd_root / "tables" / f"{profile}_neuron_representative_frames.csv"
    mapping_path = sommd_root / "inputs" / profile / "activation_core" / "activation_residue_mapping.csv"
    out_dir = sommd_root / "pymol" / profile / "neuron_representatives"
    out_dir.mkdir(parents=True, exist_ok=True)

    reps = pd.read_csv(table_path)
    mapping = pd.read_csv(mapping_path)
    reps["representative_pdb"] = ""
    reps["object_name"] = ""

    exported = 0
    selected = reps[reps["representative_som_frame"].notna()].copy()
    for (pdb_rel, traj_rel), group in selected.groupby(["pdb", "downsampled_trajectory"], sort=False):
        universe = mda.Universe(str(root / str(pdb_rel)), str(root / str(traj_rel)))
        for idx, row in group.iterrows():
            neuron = int(row["neuron"])
            cluster = int(row["cluster"])
            local_frame_index = int(row["local_frame_index"])
            universe.trajectory[local_frame_index - 1]
            object_name = f"neuron_{neuron:02d}_C{cluster}_{row['system']}_{row['replicate']}"
            out_pdb = out_dir / f"{object_name}_frame{local_frame_index:05d}.pdb"
            universe.atoms.write(str(out_pdb))
            reps.loc[idx, "representative_pdb"] = relative(out_pdb, root)
            reps.loc[idx, "object_name"] = object_name
            exported += 1

    reps.to_csv(table_path, index=False)
    pml_path = write_pymol_script(root, profile, reps, mapping)
    return table_path, pml_path, exported


def write_pymol_script(root: Path, profile: str, reps: pd.DataFrame, mapping: pd.DataFrame) -> Path:
    out_pymol = root / "article/analysis/sommd/pymol" / profile
    pml_path = out_pymol / f"{profile}_load_neuron_representatives.pml"
    lines = [
        "reinitialize",
        "bg_color white",
        "set cartoon_transparency, 0.18",
        "set sphere_scale, 0.36",
        "set label_size, 10",
    ]
    loaded: list[str] = []
    for row in reps[reps["representative_pdb"].astype(str) != ""].itertuples(index=False):
        obj = row.object_name
        system = row.system
        raw = mapping[mapping["system"] == system]["raw_resid"].astype(int).tolist()
        resid_expr = "+".join(str(x) for x in raw)
        pdb_path = root / row.representative_pdb
        lines.extend(
            [
                f"load {pdb_path.as_posix()}, {obj}",
                f"hide everything, {obj}",
                f"show cartoon, {obj}",
                f"color {SYSTEM_COLORS.get(system, 'gray50')}, {obj}",
                f"select activation_core_{obj}, {obj} and resi {resid_expr} and name CA",
                f"show spheres, activation_core_{obj}",
                f"color yellow, activation_core_{obj}",
                f"label activation_core_{obj}, resn + resi",
                f"group neuron_representatives, {obj}",
            ]
        )
        loaded.append(obj)
    if loaded:
        first = loaded[0]
        lines.extend(
            [
                "disable all",
                f"enable {first}",
                f"orient {first}",
                f"zoom {first}",
            ]
        )
    lines.append("set ray_trace_mode, 1")
    pml_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return pml_path


def main() -> int:
    args = parse_args()
    root = Path(args.root).resolve()
    profile = args.profile
    table_path = root / "article/analysis/sommd/tables" / f"{profile}_neuron_representative_frames.csv"
    if not args.skip_select:
        run_selection(root, profile, Path(args.r_prefix), table_path)
    final_table, pml_path, exported = export_pdbs(root, profile)
    print(final_table)
    print(pml_path)
    print(f"exported_pdbs={exported}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
