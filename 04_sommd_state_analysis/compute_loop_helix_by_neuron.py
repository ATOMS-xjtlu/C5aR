#!/usr/bin/env python3
from __future__ import annotations

import argparse
import warnings
from pathlib import Path

import mdtraj as md
import pandas as pd


ROOT = Path("/mnt/e/work/modeling/c5ar/md")
SOMMD_ROOT = ROOT / "article/analysis/sommd"
SYSTEM_ORDER = ("apo", "BM213", "C5apep")
SEGMENTS = {
    "ICL2": {
        "apo": range(107, 138),
        "BM213": range(95, 126),
        "C5apep": range(97, 128),
    },
    "ICL3": {
        "apo": range(204, 222),
        "BM213": range(192, 210),
        "C5apep": range(194, 212),
    },
    "H8": {
        "apo": range(283, 294),
        "BM213": range(271, 282),
        "C5apep": range(273, 284),
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Map ICL2/ICL3/H8 helix fraction onto formal SOM neurons.")
    parser.add_argument("--profile", default="formal", choices=("formal",))
    parser.add_argument("--root", default=str(ROOT))
    parser.add_argument("--chunk", type=int, default=250)
    return parser.parse_args()


def residue_indices(topology: md.Topology, system: str, segment: str) -> list[int]:
    wanted = set(SEGMENTS[segment][system])
    indices = [res.index for res in topology.residues if res.resSeq in wanted]
    found = {res.resSeq for res in topology.residues if res.resSeq in wanted}
    missing = sorted(wanted - found)
    if missing:
        raise ValueError(f"{system} {segment}: missing residue IDs {missing}")
    return indices


def compute_group(group: pd.DataFrame, root: Path, chunk: int) -> pd.DataFrame:
    first = group.iloc[0]
    system = str(first["system"])
    pdb = root / str(first["pdb"])
    traj_path = root / str(first["downsampled_trajectory"])
    by_local = group.set_index(group["local_frame_index"].astype(int))

    rows: list[pd.DataFrame] = []
    offset = 1
    for traj in md.iterload(str(traj_path), chunk=chunk, top=str(pdb)):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ss = md.compute_dssp(traj, simplified=True)
        values: dict[str, list[float]] = {}
        for segment in SEGMENTS:
            idx = residue_indices(traj.topology, system, segment)
            values[f"{segment}_helix_fraction"] = (ss[:, idx] == "H").mean(axis=1)
        local_indices = list(range(offset, offset + traj.n_frames))
        offset += traj.n_frames
        sub = by_local.reindex(local_indices)
        if sub["som_frame"].isna().any():
            missing = [local_indices[i] for i, ok in enumerate(sub["som_frame"].notna()) if not ok][:5]
            raise ValueError(f"{traj_path}: local frame indices missing from frame_mapping, first missing={missing}")
        out = sub.reset_index(drop=True).copy()
        for key, val in values.items():
            out[key] = val
        rows.append(out)

    return pd.concat(rows, ignore_index=True)


def main() -> int:
    args = parse_args()
    root = Path(args.root).resolve()
    run_dir = root / "article/analysis/sommd/runs" / args.profile / "activation_core"
    table_dir = root / "article/analysis/sommd/tables"
    table_dir.mkdir(parents=True, exist_ok=True)

    frame_map = pd.read_csv(run_dir / "frame_mapping.csv")
    neuron_clusters = pd.read_csv(run_dir / "neuron_clusters.csv")

    frame_rows: list[pd.DataFrame] = []
    group_cols = ["system", "replicate", "pdb", "downsampled_trajectory"]
    for keys, group in frame_map.groupby(group_cols, sort=False, observed=True):
        system, replicate, _, traj_path = keys
        print(f"computing DSSP: {system} {replicate} {traj_path}")
        frame_rows.append(compute_group(group.sort_values("local_frame_index"), root, args.chunk))

    frame_table = pd.concat(frame_rows, ignore_index=True)
    frame_table = frame_table.sort_values("som_frame").reset_index(drop=True)
    frame_out = table_dir / f"{args.profile}_loop_helix_by_frame.csv"
    keep_cols = [
        "global_frame",
        "profile",
        "system",
        "replicate",
        "local_frame_index",
        "time_ns",
        "som_frame",
        "neuron",
        "cluster",
    ] + [f"{segment}_helix_fraction" for segment in SEGMENTS]
    frame_table[keep_cols].to_csv(frame_out, index=False)

    neuron_rows: list[dict[str, object]] = []
    system_neuron_rows: list[dict[str, object]] = []
    system_rows: list[dict[str, object]] = []
    for segment in SEGMENTS:
        col = f"{segment}_helix_fraction"
        by_neuron = frame_table.groupby("neuron", observed=True)[col].agg(["count", "mean", "std"]).reset_index()
        by_neuron = neuron_clusters.merge(by_neuron, on="neuron", how="left")
        for row in by_neuron.itertuples(index=False):
            neuron_rows.append(
                {
                    "profile": args.profile,
                    "segment": segment,
                    "neuron": int(row.neuron),
                    "cluster": int(row.cluster),
                    "frame_count": int(row.count) if pd.notna(row.count) else 0,
                    "mean_helix_fraction": float(row.mean) if pd.notna(row.mean) else pd.NA,
                    "sd_helix_fraction": float(row.std) if pd.notna(row.std) else pd.NA,
                }
            )
        by_system_neuron = (
            frame_table.groupby(["system", "neuron"], observed=True)[col]
            .agg(["count", "mean", "std"])
            .reset_index()
            .merge(neuron_clusters, on="neuron", how="left")
        )
        for row in by_system_neuron.itertuples(index=False):
            system_neuron_rows.append(
                {
                    "profile": args.profile,
                    "segment": segment,
                    "system": row.system,
                    "neuron": int(row.neuron),
                    "cluster": int(row.cluster),
                    "frame_count": int(row.count),
                    "mean_helix_fraction": float(row.mean),
                    "sd_helix_fraction": float(row.std) if pd.notna(row.std) else pd.NA,
                }
            )
        by_system = frame_table.groupby("system", observed=True)[col].agg(["count", "mean", "std"]).reset_index()
        for row in by_system.itertuples(index=False):
            system_rows.append(
                {
                    "profile": args.profile,
                    "segment": segment,
                    "system": row.system,
                    "frame_count": int(row.count),
                    "mean_helix_fraction": float(row.mean),
                    "sd_helix_fraction": float(row.std) if pd.notna(row.std) else pd.NA,
                }
            )

    neuron_out = table_dir / f"{args.profile}_loop_helix_by_neuron.csv"
    system_neuron_out = table_dir / f"{args.profile}_loop_helix_by_system_neuron.csv"
    system_out = table_dir / f"{args.profile}_loop_helix_by_system.csv"
    pd.DataFrame(neuron_rows).to_csv(neuron_out, index=False)
    pd.DataFrame(system_neuron_rows).to_csv(system_neuron_out, index=False)
    pd.DataFrame(system_rows).to_csv(system_out, index=False)

    print(frame_out)
    print(neuron_out)
    print(system_neuron_out)
    print(system_out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
