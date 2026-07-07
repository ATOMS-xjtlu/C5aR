#!/usr/bin/env python3
"""Trim formal neuron representative PDBs to the shared 281-residue C5AR sequence."""

from __future__ import annotations

import csv
import shutil
from dataclasses import dataclass
from pathlib import Path


ROOT = Path("/mnt/e/work/modeling/c5ar/md")
FORMAL_PYMOL = ROOT / "article/analysis/sommd/pymol/formal"
NEURON_DIR = FORMAL_PYMOL / "neuron_representatives"
BACKUP_DIR = FORMAL_PYMOL / "neuron_representatives_before_sequence_unification"
TARGET_PDB = FORMAL_PYMOL / "cluster_01_BM213_rep3.pdb"
TABLE_OUT = ROOT / "article/analysis/sommd/tables/formal_neuron_representative_sequence_unification.csv"

THREE_TO_ONE = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
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
class Residue:
    chain: str
    resseq: str
    icode: str
    resname: str

    @property
    def key(self) -> tuple[str, str, str, str]:
        return self.chain, self.resseq, self.icode, self.resname


def read_residues(path: Path) -> list[Residue]:
    residues: list[Residue] = []
    seen: set[tuple[str, str, str, str]] = set()
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        residue = Residue(line[21], line[22:26].strip(), line[26].strip(), line[17:20].strip())
        if residue.key not in seen:
            seen.add(residue.key)
            residues.append(residue)
    return residues


def sequence(residues: list[Residue]) -> str:
    return "".join(THREE_TO_ONE.get(res.resname, "X") for res in residues)


def rewrite_pdb(path: Path, keep_map: dict[tuple[str, str, str, str], int]) -> None:
    lines_out: list[str] = []
    serial = 1
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        residue = Residue(line[21], line[22:26].strip(), line[26].strip(), line[17:20].strip())
        if residue.key not in keep_map:
            continue
        new_resseq = keep_map[residue.key]
        if len(line) < 80:
            line = line.ljust(80)
        new_line = (
            line[:6]
            + f"{serial:5d}"
            + line[11:21]
            + "X"
            + f"{new_resseq:4d}"
            + " "
            + line[27:]
        )
        lines_out.append(new_line.rstrip())
        serial += 1
    lines_out.append(f"TER   {serial:5d}      {lines_out[-1][17:20]} X{281:4d}")
    lines_out.append("END")
    path.write_text("\n".join(lines_out) + "\n", encoding="utf-8")


def main() -> int:
    target_residues = read_residues(TARGET_PDB)
    target_sequence = sequence(target_residues)
    if len(target_residues) != 281:
        raise ValueError(f"Expected 281 target residues, found {len(target_residues)} in {TARGET_PDB}")

    BACKUP_DIR.mkdir(parents=True, exist_ok=True)
    TABLE_OUT.parent.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []

    for pdb in sorted(NEURON_DIR.glob("neuron_*.pdb")):
        backup = BACKUP_DIR / pdb.name
        if not backup.is_file():
            shutil.copy2(pdb, backup)
        residues = read_residues(pdb)
        seq = sequence(residues)
        start = seq.find(target_sequence)
        if start < 0:
            raise ValueError(f"Target sequence not found in {pdb.name}")
        keep_residues = residues[start : start + len(target_sequence)]
        keep_map = {res.key: idx for idx, res in enumerate(keep_residues, start=1)}
        rewrite_pdb(pdb, keep_map)
        new_residues = read_residues(pdb)
        new_seq = sequence(new_residues)
        if new_seq != target_sequence:
            raise ValueError(f"Rewrite failed sequence check for {pdb.name}")
        rows.append(
            {
                "pdb": pdb.name,
                "original_residues": len(residues),
                "target_start_1based": start + 1,
                "removed_n_terminal_residues": start,
                "removed_c_terminal_residues": len(residues) - start - len(target_sequence),
                "final_residues": len(new_residues),
                "chain": "X",
                "residue_numbering": "1-281",
                "backup": str(backup),
            }
        )
        print(
            f"{pdb.name}: {len(residues)} -> {len(new_residues)} residues, "
            f"trim N={start}, C={len(residues) - start - len(target_sequence)}"
        )

    with TABLE_OUT.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {TABLE_OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
