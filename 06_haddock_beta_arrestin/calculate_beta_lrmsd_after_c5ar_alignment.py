#!/usr/bin/env python3
"""Calculate beta-arrestin CA RMSD after C5AR alignment for HADDOCK models."""

from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np


ROOT = Path("/mnt/e/work/modeling/c5ar/md")
DOCK_ROOT = ROOT / "article/analysis/docking/haddock"
TABLE_DIR = DOCK_ROOT / "tables"
MODEL_TABLE = TABLE_DIR / "haddock_cluster_model_scores.csv"
MODEL_OUT = TABLE_DIR / "haddock_beta_lrmsd_after_c5ar_alignment_model_scores.csv"
SUMMARY_OUT = TABLE_DIR / "haddock_beta_lrmsd_after_c5ar_alignment_summary.csv"

RECEPTOR_CHAIN = "X"
LIGAND_CHAIN = "A"


def read_ca(path: Path, chain: str) -> dict[int, np.ndarray]:
    coords: dict[int, np.ndarray] = {}
    with path.open() as handle:
        for line in handle:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            if line[21].strip() != chain:
                continue
            if line[12:16].strip() != "CA":
                continue
            resid_text = line[22:26].strip()
            if not resid_text.lstrip("-").isdigit():
                continue
            coords[int(resid_text)] = np.array(
                [float(line[30:38]), float(line[38:46]), float(line[46:54])],
                dtype=float,
            )
    return coords


def kabsch_fit(source: np.ndarray, target: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    source_centroid = source.mean(axis=0)
    target_centroid = target.mean(axis=0)
    source_centered = source - source_centroid
    target_centered = target - target_centroid
    covariance = source_centered.T @ target_centered
    left, _singular_values, right_t = np.linalg.svd(covariance)
    rotation = left @ right_t
    if np.linalg.det(rotation) < 0:
        left[:, -1] *= -1
        rotation = left @ right_t
    transformed = (source - source_centroid) @ rotation + target_centroid
    rmsd = float(np.sqrt(np.mean(np.sum((transformed - target) ** 2, axis=1))))
    return rotation, source_centroid, target_centroid, rmsd


def apply_transform(coords: np.ndarray, rotation: np.ndarray, source_centroid: np.ndarray, target_centroid: np.ndarray) -> np.ndarray:
    return (coords - source_centroid) @ rotation + target_centroid


def rmsd(a: np.ndarray, b: np.ndarray) -> float:
    return float(math.sqrt(float(np.mean(np.sum((a - b) ** 2, axis=1)))))


def aligned_beta_lrmsd(model_pdb: Path, reference_pdb: Path) -> dict[str, float | int | str]:
    model_receptor = read_ca(model_pdb, RECEPTOR_CHAIN)
    reference_receptor = read_ca(reference_pdb, RECEPTOR_CHAIN)
    receptor_common = sorted(set(model_receptor) & set(reference_receptor))
    if len(receptor_common) < 5:
        raise ValueError(f"Too few common receptor CA atoms: {model_pdb} vs {reference_pdb}")

    source_receptor = np.array([model_receptor[resid] for resid in receptor_common], dtype=float)
    target_receptor = np.array([reference_receptor[resid] for resid in receptor_common], dtype=float)
    rotation, source_centroid, target_centroid, receptor_fit_rmsd = kabsch_fit(source_receptor, target_receptor)

    model_ligand = read_ca(model_pdb, LIGAND_CHAIN)
    reference_ligand = read_ca(reference_pdb, LIGAND_CHAIN)
    ligand_common = sorted(set(model_ligand) & set(reference_ligand))
    if len(ligand_common) < 5:
        raise ValueError(f"Too few common beta-arrestin CA atoms: {model_pdb} vs {reference_pdb}")

    source_ligand = np.array([model_ligand[resid] for resid in ligand_common], dtype=float)
    target_ligand = np.array([reference_ligand[resid] for resid in ligand_common], dtype=float)
    fitted_ligand = apply_transform(source_ligand, rotation, source_centroid, target_centroid)
    beta_lrmsd = rmsd(fitted_ligand, target_ligand)

    return {
        "receptor_common_ca": len(receptor_common),
        "beta_common_ca": len(ligand_common),
        "receptor_ca_fit_rmsd": receptor_fit_rmsd,
        "beta_ca_lrmsd_after_c5ar_alignment": beta_lrmsd,
    }


def mean(values: list[float]) -> float:
    return sum(values) / len(values) if values else float("nan")


def stdev(values: list[float]) -> float:
    if len(values) < 2:
        return 0.0
    m = mean(values)
    return math.sqrt(sum((value - m) ** 2 for value in values) / (len(values) - 1))


def main() -> None:
    if not MODEL_TABLE.exists():
        raise FileNotFoundError(MODEL_TABLE)

    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    with MODEL_TABLE.open() as handle:
        input_rows = list(csv.DictReader(handle))

    model_rows: list[dict[str, object]] = []
    for row in input_rows:
        model_pdb = Path(row["model_pdb"])
        reference_pdb = Path(row["native_pdb"])
        metrics = aligned_beta_lrmsd(model_pdb, reference_pdb)
        model_rows.append(
            {
                "cluster": int(row["cluster"]),
                "system": row["system"],
                "model_rank": int(row["model_rank"]),
                "model_pdb": str(model_pdb),
                "haddock_score": float(row["haddock_score"]),
                "dockq": float(row["dockq"]),
                "reference_pdb": str(reference_pdb),
                "reference_mode": row.get("scoring_mode", "pseudo_native_receptor_specific"),
                "receptor_chain": RECEPTOR_CHAIN,
                "ligand_chain": LIGAND_CHAIN,
                **metrics,
            }
        )

    model_fields = [
        "cluster",
        "system",
        "model_rank",
        "model_pdb",
        "haddock_score",
        "dockq",
        "reference_pdb",
        "reference_mode",
        "receptor_chain",
        "ligand_chain",
        "receptor_common_ca",
        "beta_common_ca",
        "receptor_ca_fit_rmsd",
        "beta_ca_lrmsd_after_c5ar_alignment",
    ]
    with MODEL_OUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=model_fields)
        writer.writeheader()
        writer.writerows(model_rows)

    summary_rows: list[dict[str, object]] = []
    clusters = sorted({int(row["cluster"]) for row in model_rows})
    for cluster in clusters:
        rows = sorted((row for row in model_rows if int(row["cluster"]) == cluster), key=lambda item: int(item["model_rank"]))
        top5 = rows[:5]
        values = [float(row["beta_ca_lrmsd_after_c5ar_alignment"]) for row in top5]
        receptor_values = [float(row["receptor_ca_fit_rmsd"]) for row in top5]
        best_row = min(top5, key=lambda row: float(row["beta_ca_lrmsd_after_c5ar_alignment"]))
        summary_rows.append(
            {
                "cluster": cluster,
                "system": top5[0]["system"],
                "n_models": len(top5),
                "mean_top5_beta_lrmsd_after_c5ar_alignment": mean(values),
                "median_top5_beta_lrmsd_after_c5ar_alignment": sorted(values)[len(values) // 2],
                "sd_top5_beta_lrmsd_after_c5ar_alignment": stdev(values),
                "min_top5_beta_lrmsd_after_c5ar_alignment": min(values),
                "max_top5_beta_lrmsd_after_c5ar_alignment": max(values),
                "best_lrmsd_model_rank": best_row["model_rank"],
                "best_lrmsd_model_pdb": best_row["model_pdb"],
                "haddock_score_at_min_top5_beta_lrmsd": best_row["haddock_score"],
                "dockq_at_min_top5_beta_lrmsd": best_row["dockq"],
                "mean_top5_receptor_ca_fit_rmsd": mean(receptor_values),
                "reference_pdb": top5[0]["reference_pdb"],
                "reference_mode": top5[0]["reference_mode"],
                "receptor_chain": RECEPTOR_CHAIN,
                "ligand_chain": LIGAND_CHAIN,
            }
        )

    summary_fields = [
        "cluster",
        "system",
        "n_models",
        "mean_top5_beta_lrmsd_after_c5ar_alignment",
        "median_top5_beta_lrmsd_after_c5ar_alignment",
        "sd_top5_beta_lrmsd_after_c5ar_alignment",
        "min_top5_beta_lrmsd_after_c5ar_alignment",
        "max_top5_beta_lrmsd_after_c5ar_alignment",
        "best_lrmsd_model_rank",
        "best_lrmsd_model_pdb",
        "haddock_score_at_min_top5_beta_lrmsd",
        "dockq_at_min_top5_beta_lrmsd",
        "mean_top5_receptor_ca_fit_rmsd",
        "reference_pdb",
        "reference_mode",
        "receptor_chain",
        "ligand_chain",
    ]
    with SUMMARY_OUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_fields)
        writer.writeheader()
        writer.writerows(summary_rows)

    print(f"Wrote {MODEL_OUT}")
    print(f"Wrote {SUMMARY_OUT}")


if __name__ == "__main__":
    main()
