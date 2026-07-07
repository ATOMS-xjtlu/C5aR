#!/usr/bin/env python3
from __future__ import annotations

import argparse
import itertools
import json
from dataclasses import dataclass
from pathlib import Path

import MDAnalysis as mda
import numpy as np
import pandas as pd


ROOT = Path("/mnt/e/work/modeling/c5ar/md")
SOMMD_ROOT = ROOT / "article/analysis/sommd"
ALLOVIZ_ROOT = ROOT / "article/analysis/alloviz"

SYSTEM_ORDER = ("apo", "BM213", "C5apep")
ACTIVATION_RESIDUES = (
    "I91",
    "W102",
    "I116",
    "M120",
    "R134",
    "S171",
    "Y222",
    "F251",
    "W255",
    "Y258",
    "Y290",
    "N292",
    "N296",
    "Y300",
)
REGION_BY_LABEL = {
    "I91": "IWI / ligand pocket",
    "W102": "IWI / ligand pocket",
    "I116": "IWI / ligand pocket",
    "M120": "TM3 / DRY-like",
    "R134": "TM3 / DRY-like",
    "S171": "TM4 pocket",
    "Y222": "TM5 / DRY coupling",
    "F251": "TM6 / CWxP / PIF",
    "W255": "TM6 / CWxP / PIF",
    "Y258": "TM6 / CWxP / PIF",
    "Y290": "TM7 / NPxxY",
    "N292": "TM7 / NPxxY",
    "N296": "TM7 / NPxxY",
    "Y300": "TM7 / NPxxY",
}


@dataclass(frozen=True)
class ResidueFeature:
    feature_name: str
    residue_a: str
    residue_b: str
    canonical_resid_a: int
    canonical_resid_b: int
    region_a: str
    region_b: str
    unit: str = "angstrom"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build C5aR1 activation-core CA-CA distance features for SOMMD."
    )
    parser.add_argument(
        "--profile",
        choices=("smoke", "pilot", "formal"),
        required=True,
        help="Downsample profile to read from article/analysis/alloviz/inputs/downsampled.",
    )
    parser.add_argument(
        "--out",
        default=None,
        help="Output directory. Defaults to article/analysis/sommd/inputs/<profile>/activation_core.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing feature files.",
    )
    return parser.parse_args()


def resolve_outdir(profile: str, out: str | None) -> Path:
    if out:
        path = Path(out)
        return path if path.is_absolute() else ROOT / path
    return SOMMD_ROOT / "inputs" / profile / "activation_core"


def load_mapping() -> pd.DataFrame:
    mapping_path = ALLOVIZ_ROOT / "tables/residue_mapping_alloviz.csv"
    if not mapping_path.exists():
        raise FileNotFoundError(mapping_path)
    mapping = pd.read_csv(mapping_path)
    required = {
        "system",
        "raw_resid",
        "canonical_resid",
        "canonical_label",
        "biological_residue",
        "gpcrdb",
        "panel_context",
    }
    missing = required.difference(mapping.columns)
    if missing:
        raise ValueError(f"Residue mapping missing columns: {sorted(missing)}")
    return mapping


def canonical_ids(mapping: pd.DataFrame) -> dict[str, int]:
    ids: dict[str, int] = {}
    for label in ACTIVATION_RESIDUES:
        rows = mapping[mapping["biological_residue"].astype(str) == label]
        if rows.empty:
            rows = mapping[mapping["canonical_label"].astype(str) == label]
        found = sorted(set(rows["canonical_resid"].astype(int)))
        if len(found) != 1:
            raise ValueError(f"{label}: expected one canonical residue, found {found}")
        ids[label] = found[0]
    return ids


def residue_map(mapping: pd.DataFrame, ids: dict[str, int]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for system in SYSTEM_ORDER:
        sub = mapping[mapping["system"] == system]
        for label, canonical_resid in ids.items():
            hit = sub[sub["canonical_resid"].astype(int) == canonical_resid]
            if len(hit) != 1:
                raise ValueError(
                    f"{system} {label} canonical {canonical_resid}: expected one mapping row, found {len(hit)}"
                )
            row = hit.iloc[0]
            rows.append(
                {
                    "system": system,
                    "label": label,
                    "canonical_resid": canonical_resid,
                    "raw_resid": int(row["raw_resid"]),
                    "resname": row["resname"],
                    "gpcrdb": row.get("gpcrdb", ""),
                    "panel_context": row.get("panel_context", ""),
                    "region": REGION_BY_LABEL[label],
                }
            )
    return pd.DataFrame(rows)


def build_feature_metadata(ids: dict[str, int]) -> pd.DataFrame:
    features: list[ResidueFeature] = []
    for a, b in itertools.combinations(ACTIVATION_RESIDUES, 2):
        features.append(
            ResidueFeature(
                feature_name=f"d_{a}_{b}",
                residue_a=a,
                residue_b=b,
                canonical_resid_a=ids[a],
                canonical_resid_b=ids[b],
                region_a=REGION_BY_LABEL[a],
                region_b=REGION_BY_LABEL[b],
            )
        )
    return pd.DataFrame([f.__dict__ for f in features])


def load_manifest(profile: str) -> pd.DataFrame:
    manifest_path = ALLOVIZ_ROOT / "inputs/downsampled" / profile / "downsample_manifest.csv"
    if not manifest_path.exists():
        raise FileNotFoundError(manifest_path)
    manifest = pd.read_csv(manifest_path)
    manifest["system_order"] = manifest["system"].map({s: i for i, s in enumerate(SYSTEM_ORDER)})
    manifest["replicate_order"] = manifest["replicate"].str.extract(r"(\d+)").astype(int)
    manifest = manifest.sort_values(["system_order", "replicate_order"]).reset_index(drop=True)
    return manifest.drop(columns=["system_order", "replicate_order"])


def ca_atom_group(universe: mda.Universe, raw_resid: int, label: str, system: str):
    atoms = universe.select_atoms(f"resid {raw_resid} and name CA")
    if atoms.n_atoms != 1:
        raise ValueError(
            f"{system} {label} raw resid {raw_resid}: expected one CA atom, found {atoms.n_atoms}"
        )
    return atoms


def extract_run_features(
    row: pd.Series,
    residue_table: pd.DataFrame,
    feature_table: pd.DataFrame,
    global_start: int,
) -> tuple[np.ndarray, list[dict[str, object]]]:
    system = str(row["system"])
    replicate = str(row["replicate"])
    pdb = ROOT / str(row["pdb"])
    traj = ROOT / str(row["downsampled_trajectory"])
    if not pdb.exists():
        raise FileNotFoundError(pdb)
    if not traj.exists():
        raise FileNotFoundError(traj)

    universe = mda.Universe(str(pdb), str(traj))
    system_residues = residue_table[residue_table["system"] == system].set_index("label")
    atom_groups = [
        ca_atom_group(universe, int(system_residues.loc[label, "raw_resid"]), label, system)
        for label in ACTIVATION_RESIDUES
    ]
    pair_indices = [
        (ACTIVATION_RESIDUES.index(a), ACTIVATION_RESIDUES.index(b))
        for a, b in zip(feature_table["residue_a"], feature_table["residue_b"])
    ]

    rows: list[np.ndarray] = []
    metadata: list[dict[str, object]] = []
    stride = int(row["stride"])
    for local_idx, ts in enumerate(universe.trajectory, start=1):
        coords = np.vstack([ag.positions[0] for ag in atom_groups])
        distances = np.array(
            [np.linalg.norm(coords[i] - coords[j]) for i, j in pair_indices],
            dtype=np.float64,
        )
        rows.append(distances)
        metadata.append(
            {
                "global_frame": global_start + local_idx - 1,
                "profile": row["profile"],
                "system": system,
                "replicate": replicate,
                "pdb": row["pdb"],
                "source_trajectory": row["source_trajectory"],
                "downsampled_trajectory": row["downsampled_trajectory"],
                "stride": stride,
                "local_frame_index": local_idx,
                "source_frame_index_estimate": (local_idx - 1) * stride + 1,
                "time_ps": float(ts.time),
                "time_ns": float(ts.time) / 1000.0,
            }
        )
    matrix = np.vstack(rows)
    if not np.isfinite(matrix).all():
        raise ValueError(f"{system} {replicate}: feature matrix contains non-finite values")
    return matrix, metadata


def write_outputs(
    outdir: Path,
    manifest: pd.DataFrame,
    residue_table: pd.DataFrame,
    feature_table: pd.DataFrame,
    matrix: np.ndarray,
    metadata: pd.DataFrame,
) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    feature_names = feature_table["feature_name"].tolist()
    pd.DataFrame(matrix, columns=feature_names).to_csv(outdir / "feature_matrix.csv", index=False)
    metadata.to_csv(outdir / "frame_metadata.csv", index=False)
    feature_table.to_csv(outdir / "feature_metadata.csv", index=False)
    residue_table.to_csv(outdir / "activation_residue_mapping.csv", index=False)
    manifest.to_csv(outdir / "input_manifest_used.csv", index=False)
    summary = {
        "frames": int(matrix.shape[0]),
        "features": int(matrix.shape[1]),
        "residues": len(ACTIVATION_RESIDUES),
        "residue_labels": list(ACTIVATION_RESIDUES),
        "systems": list(SYSTEM_ORDER),
        "unit": "angstrom",
    }
    (outdir / "feature_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")


def main() -> int:
    args = parse_args()
    outdir = resolve_outdir(args.profile, args.out)
    if outdir.exists() and not args.force and (outdir / "feature_matrix.csv").exists():
        print(outdir / "feature_matrix.csv")
        return 0

    mapping = load_mapping()
    ids = canonical_ids(mapping)
    residues = residue_map(mapping, ids)
    features = build_feature_metadata(ids)
    manifest = load_manifest(args.profile)

    matrices: list[np.ndarray] = []
    metadata_rows: list[dict[str, object]] = []
    global_start = 1
    for _, row in manifest.iterrows():
        matrix, metadata = extract_run_features(row, residues, features, global_start)
        matrices.append(matrix)
        metadata_rows.extend(metadata)
        global_start += matrix.shape[0]
        print(f"[features] {row['system']} {row['replicate']}: {matrix.shape[0]} frames")

    full_matrix = np.vstack(matrices)
    metadata_df = pd.DataFrame(metadata_rows)
    if full_matrix.shape[1] != 91:
        raise ValueError(f"Expected 91 features, found {full_matrix.shape[1]}")
    if full_matrix.shape[0] != len(metadata_df):
        raise ValueError("Feature matrix and frame metadata row counts differ")
    write_outputs(outdir, manifest, residues, features, full_matrix, metadata_df)
    print(outdir)
    print(f"[features] complete: {full_matrix.shape[0]} frames x {full_matrix.shape[1]} features")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
