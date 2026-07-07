#!/usr/bin/env python3
"""Score HADDOCK cluster models with receptor-specific pseudo-native DockQ."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from Bio.Align import PairwiseAligner, substitution_matrices


ROOT = Path("/mnt/e/work/modeling/c5ar/md")
DOCK_ROOT = ROOT / "article/analysis/docking/haddock"
CASE_DIR = DOCK_ROOT / "cases"
TABLE_DIR = DOCK_ROOT / "tables"
FIG_DIR = DOCK_ROOT / "figures"
PSEUDO_NATIVE_DIR = DOCK_ROOT / "pseudo_native"
TEMPLATE_NATIVE = ROOT / "article/analysis/sommd/pymol/formal/β-GPR1.pdb"
DOCKQ_BIN = Path("/home/jiaxuan/miniconda3/bin/DockQ")
R_SCRIPT = Path("/mnt/e/app/conda-envs/R/bin/Rscript")
PLOT_SCRIPT = DOCK_ROOT / "scripts/plot_haddock_beta_arrestin_compatibility.R"

AA3_TO_1 = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "HSD": "H",
    "HSE": "H",
    "HSP": "H",
    "HID": "H",
    "HIE": "H",
    "HIP": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}


@dataclass(frozen=True)
class ModelEntry:
    cluster: int
    system: str
    rank: int
    pdb: Path
    haddock_score: float | None


@dataclass(frozen=True)
class PseudoNativeInfo:
    cluster: int
    path: Path
    aligned_ca: int
    alignment_identity: float
    receptor_ca_rmsd: float
    fit_method: str
    fit_residue_pairs: str


# GPR1 chain R residues mapped to C5AR chain X residues from
# article/analysis/sommd/pymol/formal/interface_contacts/gpr1_interface_to_c5ar_mapping.csv.
# Using the template-derived interface core avoids penalizing DockQ by a full-length
# GPR1-vs-C5AR receptor mismatch.
INTERFACE_FIT_PAIRS = [
    (68, 34),
    (70, 36),
    (76, 42),
    (135, 101),
    (138, 104),
    (140, 106),
    (141, 107),
    (142, 108),
    (143, 109),
    (151, 117),
    (226, 189),
    (238, 201),
    (239, 202),
    (244, 207),
    (247, 210),
    (304, 267),
    (308, 271),
    (309, 272),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--only-cluster", action="append", default=[], help="Cluster id(s), e.g. 1 or 1,2.")
    parser.add_argument("--max-models", type=int, default=0, help="Limit models per cluster; 0 means all water models.")
    parser.add_argument("--skip-plot", action="store_true", help="Only write tables, do not run the R plot script.")
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


def system_for_cluster(cluster: int) -> str:
    for pdb in (ROOT / "article/analysis/sommd/pymol/formal").glob(f"cluster_{cluster:02d}_*.pdb"):
        match = re.match(r"cluster_\d+_([^_]+)_.*\.pdb$", pdb.name)
        if match:
            return match.group(1)
    return ""


def parse_float(text: str) -> float | None:
    try:
        value = float(text)
    except Exception:
        return None
    return value if math.isfinite(value) else None


def score_from_line(line: str) -> float | None:
    low = line.lower()
    if "score" not in low and "haddock" not in low:
        return None
    for pattern in [
        r"haddock[-_\s]*score\s*[:=]?\s*(-?\d+(?:\.\d+)?)",
        r"score\s*[:=]\s*(-?\d+(?:\.\d+)?)",
    ]:
        match = re.search(pattern, low)
        if match:
            return parse_float(match.group(1))
    return None


def parse_haddock_score_from_pdb(pdb: Path) -> float | None:
    try:
        with pdb.open(errors="ignore") as handle:
            for line in handle:
                value = score_from_line(line)
                if value is not None:
                    return value
    except FileNotFoundError:
        return None
    return None


def parse_file_list_scores(file_list: Path) -> dict[str, float]:
    scores: dict[str, float] = {}
    if not file_list.exists():
        return scores
    for line in file_list.read_text(errors="ignore").splitlines():
        parts = line.replace('"', "").split()
        if len(parts) < 3:
            continue
        model = Path(parts[0].replace("PREVIT:", "")).name
        value = parse_float(parts[2])
        if value is not None:
            scores[model] = value
    return scores


def parse_stat_scores(water_dir: Path) -> dict[str, float]:
    scores: dict[str, float] = {}
    for stat in [
        water_dir / "analysis/structures_haddock-sorted.stat",
        water_dir / "structures_haddock-sorted.stat",
    ]:
        if not stat.exists():
            continue
        for line in stat.read_text(errors="ignore").splitlines():
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.replace('"', "").split()
            if len(parts) < 2 or not parts[0].endswith(".pdb"):
                continue
            value = parse_float(parts[1])
            if value is not None:
                scores[Path(parts[0]).name] = value
    return scores


def water_models(cluster: int, max_models: int = 0) -> list[ModelEntry]:
    case_dir = CASE_DIR / f"cluster_{cluster:02d}"
    water_dir = case_dir / "run1/structures/it1/water"
    file_list = water_dir / "file.list"
    system = system_for_cluster(cluster)
    scores = parse_file_list_scores(file_list)
    scores.update({k: v for k, v in parse_stat_scores(water_dir).items() if k not in scores})

    paths: list[Path] = []
    if file_list.exists():
        for raw in file_list.read_text(errors="ignore").splitlines():
            parts = raw.replace('"', "").split()
            if not parts:
                continue
            candidate = water_dir / Path(parts[0].replace("PREVIT:", "")).name
            if candidate.exists() and candidate.suffix.lower() == ".pdb":
                paths.append(candidate)
    if not paths and water_dir.exists():
        paths = sorted(p for p in water_dir.glob("*.pdb") if p.is_file())
    if max_models and max_models > 0:
        paths = paths[:max_models]

    out: list[ModelEntry] = []
    for idx, pdb in enumerate(paths, start=1):
        score = scores.get(pdb.name)
        if score is None:
            score = parse_haddock_score_from_pdb(pdb)
        out.append(ModelEntry(cluster=cluster, system=system, rank=idx, pdb=pdb, haddock_score=score))
    return out


def atom_coord(line: str) -> np.ndarray:
    return np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])], dtype=float)


def format_coord_line(line: str, coord: np.ndarray, serial: int | None = None, chain: str | None = None) -> str:
    line = line.rstrip("\n").ljust(80)
    if serial is not None:
        line = f"{line[:6]}{serial:5d}{line[11:]}"
    if chain is not None:
        line = f"{line[:21]}{chain}{line[22:]}"
        line = f"{line[:72]}{chain:>4}{line[76:]}"
    return f"{line[:30]}{coord[0]:8.3f}{coord[1]:8.3f}{coord[2]:8.3f}{line[54:80]}\n"


def chain_atom_lines(path: Path, chain: str) -> list[str]:
    lines: list[str] = []
    with path.open(errors="ignore") as handle:
        for line in handle:
            if line.startswith(("ATOM  ", "HETATM")) and line[21].strip() == chain:
                lines.append(line.rstrip("\n").ljust(80))
    return lines


def ca_residues(path: Path, chain: str) -> list[dict[str, object]]:
    residues: list[dict[str, object]] = []
    seen: set[tuple[int, str]] = set()
    with path.open(errors="ignore") as handle:
        for line in handle:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            if line[21].strip() != chain:
                continue
            if line[12:16] != " CA ":
                continue
            resid_text = line[22:26].strip()
            if not resid_text.lstrip("-").isdigit():
                continue
            resid = int(resid_text)
            icode = line[26].strip()
            key = (resid, icode)
            if key in seen:
                continue
            seen.add(key)
            resname = line[17:20].strip().upper()
            aa = AA3_TO_1.get(resname, "X")
            residues.append({"resid": resid, "icode": icode, "resname": resname, "aa": aa, "coord": atom_coord(line)})
    return residues


def alignment_pairs(template_res: list[dict[str, object]], target_res: list[dict[str, object]]) -> tuple[list[tuple[int, int]], float]:
    seq_template = "".join(str(r["aa"]) for r in template_res)
    seq_target = "".join(str(r["aa"]) for r in target_res)
    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -0.5
    alignment = aligner.align(seq_template, seq_target)[0]
    pairs: list[tuple[int, int]] = []
    matches = 0
    for template_block, target_block in zip(alignment.aligned[0], alignment.aligned[1]):
        t0, t1 = map(int, template_block)
        q0, q1 = map(int, target_block)
        block_len = min(t1 - t0, q1 - q0)
        for offset in range(block_len):
            ti = t0 + offset
            qi = q0 + offset
            if seq_template[ti] == "X" or seq_target[qi] == "X":
                continue
            pairs.append((ti, qi))
            if seq_template[ti] == seq_target[qi]:
                matches += 1
    identity = matches / len(pairs) if pairs else float("nan")
    return pairs, identity


def interface_fit_pairs(template_res: list[dict[str, object]], target_res: list[dict[str, object]]) -> tuple[list[tuple[int, int]], float, str]:
    template_by_resid = {int(r["resid"]): idx for idx, r in enumerate(template_res)}
    target_by_resid = {int(r["resid"]): idx for idx, r in enumerate(target_res)}
    pairs: list[tuple[int, int]] = []
    labels: list[str] = []
    matches = 0
    for gpr1_resid, c5ar_resid in INTERFACE_FIT_PAIRS:
        if gpr1_resid not in template_by_resid or c5ar_resid not in target_by_resid:
            continue
        ti = template_by_resid[gpr1_resid]
        qi = target_by_resid[c5ar_resid]
        pairs.append((ti, qi))
        labels.append(f"R:{gpr1_resid}->X:{c5ar_resid}")
        if template_res[ti]["aa"] == target_res[qi]["aa"]:
            matches += 1
    identity = matches / len(pairs) if pairs else float("nan")
    return pairs, identity, ";".join(labels)


def kabsch_fit(source: np.ndarray, target: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    source_centroid = source.mean(axis=0)
    target_centroid = target.mean(axis=0)
    source_centered = source - source_centroid
    target_centered = target - target_centroid
    covariance = source_centered.T @ target_centered
    left, _singular_values, right_t = np.linalg.svd(covariance)
    rotation = right_t.T @ left.T
    if np.linalg.det(rotation) < 0:
        right_t[-1, :] *= -1
        rotation = right_t.T @ left.T
    transformed = (source - source_centroid) @ rotation + target_centroid
    rmsd = float(np.sqrt(np.mean(np.sum((transformed - target) ** 2, axis=1))))
    return rotation, source_centroid, target_centroid, rmsd


def transform_coords(coords: np.ndarray, rotation: np.ndarray, source_centroid: np.ndarray, target_centroid: np.ndarray) -> np.ndarray:
    return (coords - source_centroid) @ rotation + target_centroid


def build_pseudo_native(cluster: int) -> PseudoNativeInfo:
    PSEUDO_NATIVE_DIR.mkdir(parents=True, exist_ok=True)
    out_path = PSEUDO_NATIVE_DIR / f"cluster_{cluster:02d}_pseudo_native.pdb"
    receptor_path = CASE_DIR / f"cluster_{cluster:02d}/receptor.pdb"
    if not receptor_path.exists():
        raise FileNotFoundError(receptor_path)

    template_res = ca_residues(TEMPLATE_NATIVE, "R")
    target_res = ca_residues(receptor_path, "X")
    pairs, identity, pair_label = interface_fit_pairs(template_res, target_res)
    fit_method = "gpr1_c5ar_getcontacts_interface_pairs"
    if len(pairs) < 5:
        pairs, identity = alignment_pairs(template_res, target_res)
        pair_label = "local_sequence_alignment"
        fit_method = "fallback_local_sequence_alignment"
    if len(pairs) < 5:
        raise RuntimeError(f"Too few GPR1-C5AR fit CA atoms for cluster {cluster}: {len(pairs)}")

    source = np.array([template_res[i]["coord"] for i, _j in pairs], dtype=float)
    target = np.array([target_res[j]["coord"] for _i, j in pairs], dtype=float)
    rotation, source_centroid, target_centroid, rmsd = kabsch_fit(source, target)

    receptor_lines = chain_atom_lines(receptor_path, "X")
    beta_lines = chain_atom_lines(TEMPLATE_NATIVE, "A")
    if not receptor_lines or not beta_lines:
        raise RuntimeError(f"Missing receptor or beta-arrestin atoms for pseudo-native cluster {cluster}")

    serial = 1
    with out_path.open("w") as out:
        out.write("REMARK receptor-specific pseudo-native for C5AR beta-arrestin DockQ\n")
        out.write(f"REMARK receptor source: {receptor_path}\n")
        out.write(f"REMARK beta template source: {TEMPLATE_NATIVE} chain A\n")
        out.write(f"REMARK fit method: {fit_method}\n")
        out.write(f"REMARK GPR1 chain R to C5AR chain X fit CA: {len(pairs)}\n")
        out.write(f"REMARK fit residue pairs: {pair_label}\n")
        out.write(f"REMARK fit-pair sequence identity: {identity:.4f}\n")
        out.write(f"REMARK receptor fit CA RMSD: {rmsd:.4f}\n")
        for line in receptor_lines:
            out.write(format_coord_line(line, atom_coord(line), serial=serial, chain="X"))
            serial += 1
        out.write("TER\n")
        for line in beta_lines:
            coord = transform_coords(atom_coord(line), rotation, source_centroid, target_centroid)
            out.write(format_coord_line(line, coord, serial=serial, chain="A"))
            serial += 1
        out.write("TER\nEND\n")
    return PseudoNativeInfo(
        cluster=cluster,
        path=out_path,
        aligned_ca=len(pairs),
        alignment_identity=identity,
        receptor_ca_rmsd=rmsd,
        fit_method=fit_method,
        fit_residue_pairs=pair_label,
    )


def extract_dockq_metrics(data: dict) -> dict[str, float | str | None]:
    best_result = data.get("best_result") or {}
    first = next(iter(best_result.values())) if isinstance(best_result, dict) and best_result else {}
    return {
        "dockq": parse_float(str(data.get("best_dockq", data.get("GlobalDockQ", first.get("DockQ", ""))))),
        "fnat": parse_float(str(first.get("fnat", first.get("Fnat", "")))) if isinstance(first, dict) else None,
        "irmsd": parse_float(str(first.get("iRMSD", first.get("irmsd", "")))) if isinstance(first, dict) else None,
        "lrmsd": parse_float(str(first.get("LRMSD", first.get("lrmsd", "")))) if isinstance(first, dict) else None,
        "best_mapping_str": data.get("best_mapping_str"),
    }


def run_dockq(model: Path, native: Path) -> dict[str, object]:
    errors: list[str] = []
    for mapping in ["XA:XA", "AX:AX"]:
        with tempfile.NamedTemporaryFile(prefix="dockq_", suffix=".json", delete=False) as tmp:
            tmp_path = Path(tmp.name)
        cmd = [str(DOCKQ_BIN), "--mapping", mapping, "--json", str(tmp_path), str(model), str(native)]
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if proc.returncode == 0 and tmp_path.exists():
            try:
                data = json.loads(tmp_path.read_text())
                metrics = extract_dockq_metrics(data)
                metrics.update(
                    {
                        "dockq_mapping_requested": mapping,
                        "dockq_returncode": proc.returncode,
                        "dockq_stdout": proc.stdout.strip()[-500:],
                        "dockq_stderr": proc.stderr.strip()[-500:],
                    }
                )
                tmp_path.unlink(missing_ok=True)
                return metrics
            except Exception as exc:  # noqa: BLE001
                errors.append(f"{mapping}: JSON parse failed: {exc}")
        else:
            errors.append(f"{mapping}: rc={proc.returncode}; stderr={proc.stderr.strip()[-240:]}")
        tmp_path.unlink(missing_ok=True)
    return {
        "dockq": None,
        "fnat": None,
        "irmsd": None,
        "lrmsd": None,
        "best_mapping_str": None,
        "dockq_mapping_requested": "",
        "dockq_returncode": 1,
        "dockq_stdout": "",
        "dockq_stderr": " | ".join(errors)[-1000:],
    }


def score_clusters(clusters: list[int], max_models: int = 0) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    rows: list[dict[str, object]] = []
    pseudo_rows: list[dict[str, object]] = []
    for cluster in clusters:
        pseudo = build_pseudo_native(cluster)
        pseudo_rows.append(
            {
                "cluster": cluster,
                "pseudo_native_pdb": str(pseudo.path),
                "template_native_pdb": str(TEMPLATE_NATIVE),
                "aligned_ca": pseudo.aligned_ca,
                "alignment_identity": pseudo.alignment_identity,
                "receptor_ca_rmsd": pseudo.receptor_ca_rmsd,
                "fit_method": pseudo.fit_method,
                "fit_residue_pairs": pseudo.fit_residue_pairs,
            }
        )
        models = water_models(cluster, max_models=max_models)
        if not models:
            rows.append(
                {
                    "cluster": cluster,
                    "system": system_for_cluster(cluster),
                    "model_rank": "",
                    "model_pdb": "",
                    "haddock_score": "",
                    "dockq": "",
                    "fnat": "",
                    "irmsd": "",
                    "lrmsd": "",
                    "scoring_mode": "pseudo_native_receptor_specific",
                    "native_pdb": str(pseudo.path),
                    "dockq_mapping_requested": "",
                    "best_mapping_str": "",
                    "dockq_returncode": "",
                    "dockq_stderr": "no water models found",
                }
            )
            continue
        for model in models:
            metrics = run_dockq(model.pdb, pseudo.path)
            row = {
                "cluster": model.cluster,
                "system": model.system,
                "model_rank": model.rank,
                "model_pdb": str(model.pdb),
                "haddock_score": model.haddock_score if model.haddock_score is not None else "",
                "dockq": metrics["dockq"] if metrics["dockq"] is not None else "",
                "fnat": metrics["fnat"] if metrics["fnat"] is not None else "",
                "irmsd": metrics["irmsd"] if metrics["irmsd"] is not None else "",
                "lrmsd": metrics["lrmsd"] if metrics["lrmsd"] is not None else "",
                "scoring_mode": "pseudo_native_receptor_specific",
                "native_pdb": str(pseudo.path),
                "pseudo_native_aligned_ca": pseudo.aligned_ca,
                "pseudo_native_alignment_identity": pseudo.alignment_identity,
                "pseudo_native_receptor_ca_rmsd": pseudo.receptor_ca_rmsd,
                "pseudo_native_fit_method": pseudo.fit_method,
                "pseudo_native_fit_residue_pairs": pseudo.fit_residue_pairs,
                "dockq_mapping_requested": metrics["dockq_mapping_requested"],
                "best_mapping_str": metrics["best_mapping_str"] or "",
                "dockq_returncode": metrics["dockq_returncode"],
                "dockq_stderr": metrics["dockq_stderr"],
            }
            rows.append(row)
            print(
                f"cluster_{cluster:02d} model {model.rank}: "
                f"pseudo-native DockQ={row['dockq']} score={row['haddock_score']}",
                flush=True,
            )
    return rows, pseudo_rows


def numeric(value: object) -> float | None:
    if value in ("", None):
        return None
    return parse_float(str(value))


def mean(values: list[float]) -> float | str:
    return sum(values) / len(values) if values else ""


def summarize(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    summaries: list[dict[str, object]] = []
    clusters = sorted({int(r["cluster"]) for r in rows})
    for cluster in clusters:
        cluster_rows = [r for r in rows if int(r["cluster"]) == cluster and numeric(r.get("dockq")) is not None]
        if not cluster_rows:
            summaries.append(
                {
                    "cluster": cluster,
                    "system": system_for_cluster(cluster),
                    "best_dockq": "",
                    "best_model_rank": "",
                    "best_model_pdb": "",
                    "haddock_score_at_best_dockq": "",
                    "mean_top5_dockq": "",
                    "mean_top5_haddock_score": "",
                    "mean_dockq": "",
                    "mean_haddock_score": "",
                    "n_models": 0,
                    "scoring_mode": "pseudo_native_receptor_specific",
                }
            )
            continue
        cluster_rows_by_rank = sorted(cluster_rows, key=lambda r: int(r["model_rank"]))
        top5 = cluster_rows_by_rank[:5]
        cluster_rows_by_dockq = sorted(
            cluster_rows,
            key=lambda r: (numeric(r.get("dockq")) or -1.0, -(numeric(r.get("model_rank")) or 0.0)),
            reverse=True,
        )
        best = cluster_rows_by_dockq[0]
        top5_dockqs = [numeric(r.get("dockq")) for r in top5 if numeric(r.get("dockq")) is not None]
        top5_scores = [numeric(r.get("haddock_score")) for r in top5 if numeric(r.get("haddock_score")) is not None]
        dockqs = [numeric(r.get("dockq")) for r in cluster_rows if numeric(r.get("dockq")) is not None]
        scores = [numeric(r.get("haddock_score")) for r in cluster_rows if numeric(r.get("haddock_score")) is not None]
        summaries.append(
            {
                "cluster": cluster,
                "system": best.get("system", system_for_cluster(cluster)),
                "best_dockq": best.get("dockq", ""),
                "best_model_rank": best.get("model_rank", ""),
                "best_model_pdb": best.get("model_pdb", ""),
                "haddock_score_at_best_dockq": best.get("haddock_score", ""),
                "mean_top5_dockq": mean([v for v in top5_dockqs if v is not None]),
                "mean_top5_haddock_score": mean([v for v in top5_scores if v is not None]),
                "mean_dockq": mean([v for v in dockqs if v is not None]),
                "mean_haddock_score": mean([v for v in scores if v is not None]),
                "n_models": len(cluster_rows),
                "scoring_mode": "pseudo_native_receptor_specific",
                "native_pdb": best.get("native_pdb", ""),
                "pseudo_native_aligned_ca": best.get("pseudo_native_aligned_ca", ""),
                "pseudo_native_alignment_identity": best.get("pseudo_native_alignment_identity", ""),
                "pseudo_native_receptor_ca_rmsd": best.get("pseudo_native_receptor_ca_rmsd", ""),
                "pseudo_native_fit_method": best.get("pseudo_native_fit_method", ""),
                "pseudo_native_fit_residue_pairs": best.get("pseudo_native_fit_residue_pairs", ""),
            }
        )
    return summaries


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError(f"No rows to write: {path}")
    fieldnames: list[str] = []
    for row in rows:
        for key in row.keys():
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def run_plots() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    proc = subprocess.run([str(R_SCRIPT), str(PLOT_SCRIPT)], cwd=ROOT, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"R plotting failed with return code {proc.returncode}")


def main() -> int:
    args = parse_args()
    if not DOCKQ_BIN.exists():
        raise FileNotFoundError(DOCKQ_BIN)
    if not TEMPLATE_NATIVE.exists():
        raise FileNotFoundError(TEMPLATE_NATIVE)
    only = selected_clusters(args.only_cluster)
    clusters = sorted(only) if only is not None else list(range(1, 10))
    rows, pseudo_rows = score_clusters(clusters, max_models=args.max_models)
    write_csv(TABLE_DIR / "haddock_pseudo_native_alignment_summary.csv", pseudo_rows)
    write_csv(TABLE_DIR / "haddock_cluster_model_scores.csv", rows)
    summaries = summarize(rows)
    write_csv(TABLE_DIR / "haddock_cluster_summary.csv", summaries)
    if not args.skip_plot and only is None:
        run_plots()
    return 0


if __name__ == "__main__":
    sys.exit(main())
