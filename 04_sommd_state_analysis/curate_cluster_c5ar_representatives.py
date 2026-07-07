#!/usr/bin/env python3
from __future__ import annotations

import argparse
import shutil
from collections import Counter
from pathlib import Path

import MDAnalysis as mda
import numpy as np
import pandas as pd


ROOT = Path("/mnt/e/work/modeling/c5ar/md")
SOMMD_ROOT = ROOT / "article/analysis/sommd"
SYSTEM_ORDER = ("apo", "BM213", "C5apep")
SYSTEM_COLORS = {"apo": "gray70", "BM213": "tv_blue", "C5apep": "orange"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Keep one dominant-system cluster representative per cluster and "
            "rewrite the PDBs as C5AR-only common-canonical-residue structures."
        )
    )
    parser.add_argument("--root", default=str(ROOT))
    parser.add_argument("--profile", default="formal", choices=("smoke", "pilot", "formal"))
    parser.add_argument("--archive-extras", action="store_true", default=True)
    parser.add_argument("--no-archive-extras", dest="archive_extras", action="store_false")
    return parser.parse_args()


def relative(path: Path, root: Path) -> str:
    return path.relative_to(root).as_posix()


def infer_offsets(residue_mapping: pd.DataFrame) -> dict[str, int]:
    offsets: dict[str, int] = {}
    for system, sub in residue_mapping.groupby("system", sort=False):
        values = (sub["raw_resid"].astype(int) - sub["canonical_resid"].astype(int)).tolist()
        counts = Counter(values)
        offset, n = counts.most_common(1)[0]
        if n != len(values):
            raise ValueError(f"{system}: raw-canonical residue offset is not constant: {counts}")
        offsets[str(system)] = int(offset)
    return offsets


def protein_residue_map(pdb_path: Path, offset: int) -> dict[int, tuple[int, str]]:
    universe = mda.Universe(str(pdb_path))
    protein = universe.select_atoms("protein")
    out: dict[int, tuple[int, str]] = {}
    for residue in protein.residues:
        raw_resid = int(residue.resid)
        canonical_resid = raw_resid - offset
        out[canonical_resid] = (raw_resid, str(residue.resname))
    return out


def common_c5ar_residues(
    root: Path,
    manifest: pd.DataFrame,
    offsets: dict[str, int],
) -> tuple[list[int], dict[str, dict[int, tuple[int, str]]], pd.DataFrame]:
    maps: dict[str, dict[int, tuple[int, str]]] = {}
    for system in SYSTEM_ORDER:
        pdb_rel = manifest.loc[manifest["system"] == system, "pdb"].iloc[0]
        maps[system] = protein_residue_map(root / str(pdb_rel), offsets[system])

    common = sorted(set.intersection(*(set(maps[system]) for system in SYSTEM_ORDER)))
    rows: list[dict[str, object]] = []
    for canonical_resid in common:
        row: dict[str, object] = {"canonical_resid": canonical_resid}
        names = []
        for system in SYSTEM_ORDER:
            raw_resid, resname = maps[system][canonical_resid]
            row[f"{system}_raw_resid"] = raw_resid
            row[f"{system}_resname"] = resname
            names.append(resname)
        row["resname_consistent"] = len(set(names)) == 1
        rows.append(row)
    return common, maps, pd.DataFrame(rows)


def dominant_cluster_representatives(
    frame_map: pd.DataFrame,
    representative_frames: pd.DataFrame,
) -> pd.DataFrame:
    counts = (
        frame_map.groupby(["cluster", "system"], observed=True)
        .size()
        .rename("cluster_system_frames")
        .reset_index()
    )
    cluster_totals = (
        frame_map.groupby("cluster", observed=True)
        .size()
        .rename("cluster_total_frames")
        .reset_index()
    )
    counts = counts.merge(cluster_totals, on="cluster", how="left")
    counts["fraction_within_cluster"] = counts["cluster_system_frames"] / counts["cluster_total_frames"]
    counts["system_rank"] = counts.groupby("cluster")["cluster_system_frames"].rank(method="first", ascending=False)
    dominant = counts[counts["system_rank"] == 1].drop(columns=["system_rank"])

    reps = representative_frames.merge(
        dominant,
        on=["cluster", "system"],
        how="left",
        validate="one_to_one",
    )
    missing = reps[reps["cluster_system_frames"].isna()]
    if not missing.empty:
        bad = ", ".join(f"C{int(row.cluster)}:{row.system}" for row in missing.itertuples(index=False))
        raise ValueError(f"Representative rows are not dominant-system rows: {bad}")
    reps = reps.sort_values("cluster").reset_index(drop=True)
    return reps


def rewrite_pdb_resids(tmp_pdb: Path, out_pdb: Path, raw_to_canonical: dict[int, int]) -> None:
    lines: list[str] = []
    for line in tmp_pdb.read_text(encoding="utf-8").splitlines():
        if line.startswith(("ATOM", "HETATM", "TER")) and len(line) >= 26:
            try:
                raw_resid = int(line[22:26])
            except ValueError:
                lines.append(line)
                continue
            if raw_resid not in raw_to_canonical:
                continue
            canonical_resid = raw_to_canonical[raw_resid]
            line = f"{line[:21]}X{canonical_resid:4d}{line[26:]}"
        lines.append(line)
    if not any(line.startswith("END") for line in lines):
        lines.append("END")
    out_pdb.write_text("\n".join(lines) + "\n", encoding="utf-8")


def count_pdb_atoms_residues(path: Path) -> tuple[int, int]:
    atoms = 0
    residues: set[tuple[str, int]] = set()
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(("ATOM", "HETATM")):
            atoms += 1
            chain = line[21].strip() or "_"
            resid = int(line[22:26])
            residues.add((chain, resid))
    return atoms, len(residues)


def export_c5ar_only_representatives(
    root: Path,
    representatives: pd.DataFrame,
    common_resids: list[int],
    offsets: dict[str, int],
    out_pymol: Path,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    out_pymol.mkdir(parents=True, exist_ok=True)

    for row in representatives.itertuples(index=False):
        system = str(row.system)
        offset = offsets[system]
        raw_resids = [canonical_resid + offset for canonical_resid in common_resids]
        raw_set = set(raw_resids)
        raw_to_canonical = {raw: canonical for canonical, raw in zip(common_resids, raw_resids)}

        universe = mda.Universe(str(root / str(row.pdb)), str(root / str(row.downsampled_trajectory)))
        universe.trajectory[int(row.local_frame_index) - 1]
        protein = universe.select_atoms("protein")
        mask = np.isin(protein.resids.astype(int), list(raw_set))
        selected = protein[mask]

        object_name = f"cluster_{int(row.cluster):02d}_{system}_{row.replicate}"
        out_pdb = out_pymol / f"{object_name}.pdb"
        tmp_pdb = out_pymol / f".{object_name}.tmp.pdb"
        selected.write(str(tmp_pdb))
        rewrite_pdb_resids(tmp_pdb, out_pdb, raw_to_canonical)
        tmp_pdb.unlink(missing_ok=True)
        atoms, residues = count_pdb_atoms_residues(out_pdb)

        record = row._asdict()
        record.update(
            {
                "object_name": object_name,
                "representative_pdb": relative(out_pdb, root),
                "c5ar_common_canonical_residues": len(common_resids),
                "pdb_atoms": atoms,
                "pdb_residues": residues,
                "raw_residue_min": min(raw_resids),
                "raw_residue_max": max(raw_resids),
                "canonical_residue_min": min(common_resids),
                "canonical_residue_max": max(common_resids),
            }
        )
        rows.append(record)

    return pd.DataFrame(rows)


def archive_extra_cluster_pdbs(out_pymol: Path, keep_names: set[str]) -> list[Path]:
    archive_dir = out_pymol / "archive_non_dominant_cluster_pdbs"
    moved: list[Path] = []
    for path in sorted(out_pymol.glob("cluster_*.pdb")):
        if path.name in keep_names:
            continue
        archive_dir.mkdir(parents=True, exist_ok=True)
        target = archive_dir / path.name
        if target.exists():
            stem = target.stem
            suffix = target.suffix
            i = 1
            while True:
                candidate = archive_dir / f"{stem}.bak{i}{suffix}"
                if not candidate.exists():
                    target = candidate
                    break
                i += 1
        shutil.move(str(path), str(target))
        moved.append(target)
    return moved


def write_pymol_script(
    representatives: pd.DataFrame,
    activation_mapping: pd.DataFrame,
    root: Path,
    out_pymol: Path,
    profile: str,
) -> Path:
    canonical = (
        activation_mapping[["canonical_resid"]]
        .drop_duplicates()
        .sort_values("canonical_resid")["canonical_resid"]
        .astype(int)
        .tolist()
    )
    resid_expr = "+".join(str(x) for x in canonical)
    lines = [
        "reinitialize",
        "bg_color white",
        "set cartoon_transparency, 0.15",
        "set sphere_scale, 0.45",
        "set label_size, 14",
    ]
    for row in representatives.itertuples(index=False):
        obj = row.object_name
        pdb_path = root / str(row.representative_pdb)
        lines.extend(
            [
                f"load {pdb_path.as_posix()}, {obj}",
                f"hide everything, {obj}",
                f"show cartoon, {obj}",
                f"color {SYSTEM_COLORS.get(row.system, 'gray50')}, {obj}",
                f"select activation_core_{obj}, {obj} and chain X and resi {resid_expr} and name CA",
                f"show spheres, activation_core_{obj}",
                f"color yellow, activation_core_{obj}",
                f"label activation_core_{obj}, resn + resi",
            ]
        )
    if not representatives.empty:
        first = representatives.iloc[0]["object_name"]
        for obj in representatives["object_name"].iloc[1:]:
            lines.append(f"align {obj} and chain X and name CA, {first} and chain X and name CA")
    lines.extend(["orient", "zoom", "set ray_trace_mode, 1"])
    pml = out_pymol / f"{profile}_load_cluster_representatives.pml"
    pml.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return pml


def main() -> int:
    args = parse_args()
    root = Path(args.root).resolve()
    profile = args.profile
    sommd_root = root / "article/analysis/sommd"
    run_dir = sommd_root / "runs" / profile / "activation_core"
    input_dir = sommd_root / "inputs" / profile / "activation_core"
    table_dir = sommd_root / "tables"
    out_pymol = sommd_root / "pymol" / profile
    table_dir.mkdir(parents=True, exist_ok=True)
    out_pymol.mkdir(parents=True, exist_ok=True)

    frame_map = pd.read_csv(run_dir / "frame_mapping.csv")
    representatives = pd.read_csv(table_dir / f"{profile}_representative_frames_full.csv")
    activation_mapping = pd.read_csv(input_dir / "activation_residue_mapping.csv")
    manifest = pd.read_csv(input_dir / "input_manifest_used.csv")

    offsets = infer_offsets(activation_mapping)
    common_resids, _, common_table = common_c5ar_residues(root, manifest, offsets)
    common_table.to_csv(table_dir / f"{profile}_c5ar_common_canonical_residues.csv", index=False)

    curated = dominant_cluster_representatives(frame_map, representatives)
    exported = export_c5ar_only_representatives(root, curated, common_resids, offsets, out_pymol)
    keep_names = {Path(path).name for path in exported["representative_pdb"]}
    moved = archive_extra_cluster_pdbs(out_pymol, keep_names) if args.archive_extras else []
    pml = write_pymol_script(exported, activation_mapping, root, out_pymol, profile)

    out_csv = table_dir / f"{profile}_cluster_representatives_c5ar_only.csv"
    exported.to_csv(out_csv, index=False)

    print(out_csv)
    print(pml)
    print(f"clusters={len(exported)}")
    print(f"common_c5ar_residues={len(common_resids)}")
    print(f"residue_counts={sorted(exported['pdb_residues'].unique().tolist())}")
    print(f"atom_counts={sorted(exported['pdb_atoms'].unique().tolist())}")
    print(f"archived_extra_pdbs={len(moved)}")
    for path in moved:
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
