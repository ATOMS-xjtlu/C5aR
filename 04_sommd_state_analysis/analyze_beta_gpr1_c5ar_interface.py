#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Summarize beta-arrestin/GPR1 contacts and map GPR1 interface sites to C5AR."""

from __future__ import annotations

import csv
import math
from collections import Counter, defaultdict
from pathlib import Path

from Bio import Align
from Bio.Align import substitution_matrices


ROOT = Path(__file__).resolve().parents[4]
FORMAL = ROOT / "article/analysis/sommd/pymol/formal"
OUTDIR = FORMAL / "interface_contacts"

BETA_GPR1_PDB = FORMAL / "β-GPR1.pdb"
GETCONTACTS_TSV = OUTDIR / "beta_gpr1_getcontacts_chainR_chainA.tsv"
C5AR_REPS = {
    "cluster_01_BM213_rep3": FORMAL / "cluster_01_BM213_rep3.pdb",
    "cluster_02_C5apep_rep3": FORMAL / "cluster_02_C5apep_rep3.pdb",
}
C5AR_PAPER_NUMBER_OFFSET = 33
HADDOCK_ACTIVE_C5AR_STRUCT = [108, 109, 207, 210, 271, 272]
HADDOCK_PASSIVE_C5AR_STRUCT = [36, 101, 106, 201]
HADDOCK_ACTIVE_BETA = [67, 71, 73, 244, 249]
HADDOCK_PASSIVE_BETA = [66, 69, 70, 241, 245]

AA3 = {
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
    "MSE": "M",
}


def parse_pdb_residues(pdb: Path, chain: str) -> list[dict[str, object]]:
    residues: list[dict[str, object]] = []
    seen: set[tuple[str, str, str]] = set()
    with pdb.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith("ATOM"):
                continue
            if line[21].strip() != chain:
                continue
            key = (line[21], line[22:26], line[26])
            if key in seen:
                continue
            seen.add(key)
            resname = line[17:20].strip()
            residues.append(
                {
                    "chain": chain,
                    "resname": resname,
                    "resid": int(line[22:26]),
                    "icode": line[26].strip(),
                    "aa": AA3.get(resname, "X"),
                }
            )
    return residues


def parse_pdb_atoms(pdb: Path, chains: set[str]) -> list[dict[str, object]]:
    atoms: list[dict[str, object]] = []
    with pdb.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith("ATOM"):
                continue
            chain = line[21].strip()
            if chain not in chains:
                continue
            element = (line[76:78].strip() or line[12:16].strip()[0]).upper()
            if element.startswith("H"):
                continue
            atoms.append(
                {
                    "chain": chain,
                    "resname": line[17:20].strip(),
                    "resid": int(line[22:26]),
                    "atom": line[12:16].strip(),
                    "x": float(line[30:38]),
                    "y": float(line[38:46]),
                    "z": float(line[46:54]),
                    "line": line,
                }
            )
    return atoms


def parse_getcontacts(tsv: Path):
    residue_stats = defaultdict(
        lambda: {"count": 0, "types": Counter(), "min_distance": math.inf, "partners": set()}
    )
    pair_stats = defaultdict(lambda: {"count": 0, "types": Counter(), "min_distance": math.inf})

    with tsv.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            itype = parts[1]
            atom_tokens = [part for part in parts[2:-1] if part.count(":") >= 3]
            dist = float(parts[-1])
            parsed = []
            for token in atom_tokens:
                chain, resname, resid, atom = token.split(":")[:4]
                parsed.append((chain, resname, int(resid), atom))
            receptor_atoms = [atom for atom in parsed if atom[0] == "R"]
            beta_atoms = [atom for atom in parsed if atom[0] == "A"]
            if not receptor_atoms or not beta_atoms:
                continue
            receptor = receptor_atoms[0]
            beta = beta_atoms[0]
            pkey = (receptor[0], receptor[1], receptor[2], beta[0], beta[1], beta[2])
            pair_stats[pkey]["count"] += 1
            pair_stats[pkey]["types"][itype] += 1
            pair_stats[pkey]["min_distance"] = min(pair_stats[pkey]["min_distance"], dist)
            for side, residue, partner in (
                ("gpr1", receptor, beta),
                ("beta_arrestin", beta, receptor),
            ):
                rkey = (side, residue[0], residue[1], residue[2])
                residue_stats[rkey]["count"] += 1
                residue_stats[rkey]["types"][itype] += 1
                residue_stats[rkey]["min_distance"] = min(
                    residue_stats[rkey]["min_distance"], dist
                )
                residue_stats[rkey]["partners"].add((partner[0], partner[1], partner[2]))
    return residue_stats, pair_stats


def is_recommended(stats: dict[str, object]) -> bool:
    types = stats["types"]
    return bool(set(types) - {"vdw"}) or stats["count"] >= 3 or stats["min_distance"] <= 3.5


def type_string(types: Counter) -> str:
    return ";".join(f"{key}:{types[key]}" for key in sorted(types))


def partner_string(partners: set[tuple[str, str, int]]) -> str:
    return ";".join(f"{chain}:{resname}{resid}" for chain, resname, resid in sorted(partners, key=lambda x: (x[0], x[2], x[1])))


def heavy_atom_contacts(pdb: Path, cutoff: float = 4.5):
    atoms = parse_pdb_atoms(pdb, {"R", "A"})
    receptor_atoms = [atom for atom in atoms if atom["chain"] == "R"]
    beta_atoms = [atom for atom in atoms if atom["chain"] == "A"]
    cutoff2 = cutoff * cutoff
    pairs = defaultdict(lambda: {"min_distance": math.inf, "atom_contacts": 0})
    residues = {
        "gpr1": defaultdict(lambda: {"min_distance": math.inf, "atom_contacts": 0, "partners": set()}),
        "beta_arrestin": defaultdict(
            lambda: {"min_distance": math.inf, "atom_contacts": 0, "partners": set()}
        ),
    }
    for receptor in receptor_atoms:
        rx, ry, rz = receptor["x"], receptor["y"], receptor["z"]
        for beta in beta_atoms:
            dx = rx - beta["x"]
            dy = ry - beta["y"]
            dz = rz - beta["z"]
            dist2 = dx * dx + dy * dy + dz * dz
            if dist2 > cutoff2:
                continue
            dist = math.sqrt(dist2)
            pkey = (
                receptor["chain"],
                receptor["resname"],
                receptor["resid"],
                beta["chain"],
                beta["resname"],
                beta["resid"],
            )
            pairs[pkey]["min_distance"] = min(pairs[pkey]["min_distance"], dist)
            pairs[pkey]["atom_contacts"] += 1
            for side, atom, partner in (
                ("gpr1", receptor, beta),
                ("beta_arrestin", beta, receptor),
            ):
                rkey = (atom["chain"], atom["resname"], atom["resid"])
                residues[side][rkey]["min_distance"] = min(
                    residues[side][rkey]["min_distance"], dist
                )
                residues[side][rkey]["atom_contacts"] += 1
                residues[side][rkey]["partners"].add(
                    (partner["chain"], partner["resname"], partner["resid"])
                )
    return residues, pairs


def build_sequence_mapping(gpr1_residues: list[dict[str, object]], c5ar_residues: list[dict[str, object]]):
    seq1 = "".join(str(residue["aa"]) for residue in gpr1_residues)
    seq2 = "".join(str(residue["aa"]) for residue in c5ar_residues)
    aligner = Align.PairwiseAligner(mode="global")
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    alignment = aligner.align(seq1, seq2)[0]

    mapping: dict[int, dict[str, object]] = {}
    aligned_cols = 0
    identical_cols = 0
    for block1, block2 in zip(alignment.aligned[0], alignment.aligned[1]):
        for i, j in zip(range(block1[0], block1[1]), range(block2[0], block2[1])):
            mapping[int(gpr1_residues[i]["resid"])] = c5ar_residues[j]
            aligned_cols += 1
            if seq1[i] == seq2[j]:
                identical_cols += 1
    identity = identical_cols / aligned_cols if aligned_cols else 0
    return mapping, alignment.score, aligned_cols, identity


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_site_file(path: Path, chain: str, residues: list[int]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for resid in residues:
            handle.write(f"{chain}:{resid}\n")


def one_letter(resname: str) -> str:
    return AA3.get(resname, "X")


def c5ar_paper_label(resname: str, structure_resid: int) -> str:
    return f"{one_letter(resname)}{structure_resid + C5AR_PAPER_NUMBER_OFFSET}"


def residue_lookup(residues: list[dict[str, object]]) -> dict[int, dict[str, object]]:
    return {int(row["resid"]): row for row in residues}


def evidence_for_c5ar_struct_resid(
    mapping_rows: list[dict[str, object]], rep_name: str, struct_resid: int
) -> dict[str, object] | None:
    key = f"{rep_name}_c5ar_resid"
    for row in mapping_rows:
        if row.get(key) != "" and int(row[key]) == struct_resid:
            return row
    return None


def evidence_for_beta_resid(
    residue_rows: list[dict[str, object]], beta_resid: int
) -> dict[str, object] | None:
    for row in residue_rows:
        if row["side"] == "beta_arrestin" and int(row["resid"]) == beta_resid:
            return row
    return None


def write_beta_chain_pdb(src: Path, dest: Path) -> None:
    with src.open(encoding="utf-8", errors="replace") as inp, dest.open("w", encoding="utf-8") as out:
        serial = 1
        for line in inp:
            if not line.startswith("ATOM") or line[21].strip() != "A":
                continue
            out.write(f"{line[:6]}{serial:5d}{line[11:]}")
            serial += 1
        out.write("TER\nEND\n")


def residue_sort_key(row: dict[str, object]):
    return (str(row.get("chain", "")), int(row["resid"]))


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    residue_stats, pair_stats = parse_getcontacts(GETCONTACTS_TSV)
    heavy_residues, heavy_pairs = heavy_atom_contacts(BETA_GPR1_PDB, cutoff=4.5)

    residue_rows = []
    for (side, chain, resname, resid), stats in sorted(
        residue_stats.items(), key=lambda item: (item[0][0], item[0][3])
    ):
        residue_rows.append(
            {
                "side": side,
                "chain": chain,
                "resname": resname,
                "resid": resid,
                "atom_contact_count_getcontacts": stats["count"],
                "contact_types": type_string(stats["types"]),
                "min_distance_A_getcontacts": f"{stats['min_distance']:.3f}",
                "partner_residue_count": len(stats["partners"]),
                "partner_residues": partner_string(stats["partners"]),
                "recommended_active": "yes" if is_recommended(stats) else "no",
            }
        )
    write_csv(
        OUTDIR / "beta_gpr1_getcontacts_interface_residues.csv",
        residue_rows,
        [
            "side",
            "chain",
            "resname",
            "resid",
            "atom_contact_count_getcontacts",
            "contact_types",
            "min_distance_A_getcontacts",
            "partner_residue_count",
            "partner_residues",
            "recommended_active",
        ],
    )

    pair_rows = []
    for (r_chain, r_resname, r_resid, b_chain, b_resname, b_resid), stats in sorted(
        pair_stats.items(), key=lambda item: (item[0][2], item[0][5])
    ):
        pair_rows.append(
            {
                "gpr1_chain": r_chain,
                "gpr1_resname": r_resname,
                "gpr1_resid": r_resid,
                "beta_chain": b_chain,
                "beta_resname": b_resname,
                "beta_resid": b_resid,
                "atom_contact_count_getcontacts": stats["count"],
                "contact_types": type_string(stats["types"]),
                "min_distance_A_getcontacts": f"{stats['min_distance']:.3f}",
                "recommended_pair": "yes" if is_recommended(stats) else "no",
            }
        )
    write_csv(
        OUTDIR / "beta_gpr1_getcontacts_interface_pairs.csv",
        pair_rows,
        [
            "gpr1_chain",
            "gpr1_resname",
            "gpr1_resid",
            "beta_chain",
            "beta_resname",
            "beta_resid",
            "atom_contact_count_getcontacts",
            "contact_types",
            "min_distance_A_getcontacts",
            "recommended_pair",
        ],
    )

    heavy_rows = []
    for (r_chain, r_resname, r_resid, b_chain, b_resname, b_resid), stats in sorted(
        heavy_pairs.items(), key=lambda item: (item[0][2], item[0][5])
    ):
        heavy_rows.append(
            {
                "gpr1_chain": r_chain,
                "gpr1_resname": r_resname,
                "gpr1_resid": r_resid,
                "beta_chain": b_chain,
                "beta_resname": b_resname,
                "beta_resid": b_resid,
                "heavy_atom_contact_count_4p5A": stats["atom_contacts"],
                "min_distance_A_4p5": f"{stats['min_distance']:.3f}",
            }
        )
    write_csv(
        OUTDIR / "beta_gpr1_heavy_atom_pairs_4p5A.csv",
        heavy_rows,
        [
            "gpr1_chain",
            "gpr1_resname",
            "gpr1_resid",
            "beta_chain",
            "beta_resname",
            "beta_resid",
            "heavy_atom_contact_count_4p5A",
            "min_distance_A_4p5",
        ],
    )

    gpr1_seq = parse_pdb_residues(BETA_GPR1_PDB, "R")
    mappings = {}
    alignment_notes = {}
    c5ar_residue_lookups = {}
    for name, pdb in C5AR_REPS.items():
        c5ar_seq = parse_pdb_residues(pdb, "X")
        c5ar_residue_lookups[name] = residue_lookup(c5ar_seq)
        mapping, score, aligned_cols, identity = build_sequence_mapping(gpr1_seq, c5ar_seq)
        mappings[name] = mapping
        alignment_notes[name] = {
            "score": score,
            "aligned_cols": aligned_cols,
            "identity": identity,
            "residue_count": len(c5ar_seq),
        }

    getcontacts_gpr1_resids = {
        resid
        for (side, _chain, _resname, resid), _stats in residue_stats.items()
        if side == "gpr1"
    }
    heavy_gpr1_resids = {resid for (_chain, _resname, resid) in heavy_residues["gpr1"]}
    mapping_rows = []
    for resid in sorted(getcontacts_gpr1_resids | heavy_gpr1_resids):
        gpr1_residue = next(residue for residue in gpr1_seq if residue["resid"] == resid)
        contact_key = ("gpr1", "R", gpr1_residue["resname"], resid)
        contact_stats = residue_stats.get(contact_key)
        heavy_key = ("R", gpr1_residue["resname"], resid)
        heavy_stats = heavy_residues["gpr1"].get(heavy_key)
        row = {
            "gpr1_chain": "R",
            "gpr1_resname": gpr1_residue["resname"],
            "gpr1_resid": resid,
            "in_getcontacts": "yes" if contact_stats else "no",
            "getcontacts_recommended_active": "yes" if contact_stats and is_recommended(contact_stats) else "no",
            "getcontacts_contact_types": type_string(contact_stats["types"]) if contact_stats else "",
            "getcontacts_atom_contact_count": contact_stats["count"] if contact_stats else 0,
            "getcontacts_min_distance_A": f"{contact_stats['min_distance']:.3f}" if contact_stats else "",
            "in_heavy_atom_4p5A": "yes" if heavy_stats else "no",
            "heavy_atom_contact_count_4p5A": heavy_stats["atom_contacts"] if heavy_stats else 0,
            "heavy_atom_min_distance_A_4p5": f"{heavy_stats['min_distance']:.3f}" if heavy_stats else "",
        }
        for name, mapping in mappings.items():
            mapped = mapping.get(resid)
            prefix = name
            row[f"{prefix}_c5ar_chain"] = mapped["chain"] if mapped else ""
            row[f"{prefix}_c5ar_resname"] = mapped["resname"] if mapped else ""
            row[f"{prefix}_c5ar_resid"] = mapped["resid"] if mapped else ""
        mapping_rows.append(row)

    mapping_fields = [
        "gpr1_chain",
        "gpr1_resname",
        "gpr1_resid",
        "in_getcontacts",
        "getcontacts_recommended_active",
        "getcontacts_contact_types",
        "getcontacts_atom_contact_count",
        "getcontacts_min_distance_A",
        "in_heavy_atom_4p5A",
        "heavy_atom_contact_count_4p5A",
        "heavy_atom_min_distance_A_4p5",
    ]
    for name in C5AR_REPS:
        mapping_fields += [
            f"{name}_c5ar_chain",
            f"{name}_c5ar_resname",
            f"{name}_c5ar_resid",
        ]
    write_csv(OUTDIR / "gpr1_interface_to_c5ar_mapping.csv", mapping_rows, mapping_fields)

    core_c5ar_by_rep = {}
    all_c5ar_by_rep = {}
    passive_c5ar_by_rep = {}
    for name in C5AR_REPS:
        core_c5ar_by_rep[name] = sorted(
            {
                int(row[f"{name}_c5ar_resid"])
                for row in mapping_rows
                if row["getcontacts_recommended_active"] == "yes" and row[f"{name}_c5ar_resid"] != ""
            }
        )
        all_c5ar_by_rep[name] = sorted(
            {
                int(row[f"{name}_c5ar_resid"])
                for row in mapping_rows
                if row["in_getcontacts"] == "yes" and row[f"{name}_c5ar_resid"] != ""
            }
        )
        passive_c5ar_by_rep[name] = sorted(
            {
                int(row[f"{name}_c5ar_resid"])
                for row in mapping_rows
                if row["in_heavy_atom_4p5A"] == "yes" and row[f"{name}_c5ar_resid"] != ""
            }
        )
        write_site_file(OUTDIR / f"hdock_rsite_{name}_core.txt", "X", core_c5ar_by_rep[name])
        write_site_file(OUTDIR / f"hdock_rsite_{name}_all_getcontacts.txt", "X", all_c5ar_by_rep[name])
        write_site_file(OUTDIR / f"haddock_passive_rsite_{name}_heavy4p5A.txt", "X", passive_c5ar_by_rep[name])

    beta_core = sorted(
        resid
        for (side, _chain, _resname, resid), stats in residue_stats.items()
        if side == "beta_arrestin" and is_recommended(stats)
    )
    beta_all = sorted(
        resid
        for (side, _chain, _resname, resid), _stats in residue_stats.items()
        if side == "beta_arrestin"
    )
    beta_heavy = sorted(resid for (_chain, _resname, resid) in heavy_residues["beta_arrestin"])
    write_site_file(OUTDIR / "hdock_lsite_beta_arrestin_core.txt", "A", beta_core)
    write_site_file(OUTDIR / "hdock_lsite_beta_arrestin_all_getcontacts.txt", "A", beta_all)
    write_site_file(OUTDIR / "haddock_passive_lsite_beta_arrestin_heavy4p5A.txt", "A", beta_heavy)
    for name in C5AR_REPS:
        write_site_file(OUTDIR / f"haddock_active_rsite_{name}.txt", "X", HADDOCK_ACTIVE_C5AR_STRUCT)
        write_site_file(
            OUTDIR / f"haddock_passive_rsite_{name}_conservative.txt",
            "X",
            HADDOCK_PASSIVE_C5AR_STRUCT,
        )
    write_site_file(OUTDIR / "haddock_active_lsite_beta_arrestin.txt", "A", HADDOCK_ACTIVE_BETA)
    write_site_file(
        OUTDIR / "haddock_passive_lsite_beta_arrestin_conservative.txt",
        "A",
        HADDOCK_PASSIVE_BETA,
    )

    haddock_numbering_rows = build_haddock_numbering_rows(
        residue_rows=residue_rows,
        mapping_rows=mapping_rows,
        c5ar_residue_lookups=c5ar_residue_lookups,
    )
    write_csv(
        OUTDIR / "haddock_active_passive_numbering.csv",
        haddock_numbering_rows,
        [
            "molecule",
            "role",
            "docking_chain",
            "docking_structure_resid",
            "docking_structure_residue",
            "paper_or_crystal_residue",
            "paper_or_crystal_resid",
            "source_template_residue",
            "evidence",
            "notes",
        ],
    )
    write_beta_chain_pdb(BETA_GPR1_PDB, FORMAL / "beta_arrestin_chainA_from_beta_GPR1.pdb")

    report = build_report(
        residue_rows=residue_rows,
        pair_rows=pair_rows,
        mapping_rows=mapping_rows,
        alignment_notes=alignment_notes,
        core_c5ar_by_rep=core_c5ar_by_rep,
        all_c5ar_by_rep=all_c5ar_by_rep,
        passive_c5ar_by_rep=passive_c5ar_by_rep,
        beta_core=beta_core,
        beta_all=beta_all,
        beta_heavy=beta_heavy,
        haddock_numbering_rows=haddock_numbering_rows,
    )
    (FORMAL / "beta_gpr1_c5ar_interface_constraints.md").write_text(report, encoding="utf-8")


def format_site_list(chain: str, residues: list[int]) -> str:
    return ", ".join(f"{chain}:{resid}" for resid in residues)


def format_residue_list(rows: list[dict[str, object]], side: str, recommended: bool | None = None) -> str:
    selected = [row for row in rows if row["side"] == side]
    if recommended is not None:
        selected = [row for row in selected if (row["recommended_active"] == "yes") == recommended]
    selected = sorted(selected, key=lambda row: int(row["resid"]))
    return ", ".join(f"{row['chain']}:{row['resname']}{row['resid']}" for row in selected)


def build_haddock_numbering_rows(*, residue_rows, mapping_rows, c5ar_residue_lookups):
    rows = []
    rep_name = "cluster_02_C5apep_rep3"
    lookup = c5ar_residue_lookups[rep_name]
    role_map = [
        ("active", HADDOCK_ACTIVE_C5AR_STRUCT),
        ("passive", HADDOCK_PASSIVE_C5AR_STRUCT),
    ]
    for role, struct_resids in role_map:
        for struct_resid in struct_resids:
            residue = lookup[struct_resid]
            evidence = evidence_for_c5ar_struct_resid(mapping_rows, rep_name, struct_resid)
            if evidence:
                template_residue = f"GPR1 R:{evidence['gpr1_resname']}{evidence['gpr1_resid']}"
                ev_text = (
                    f"{evidence['getcontacts_contact_types']}; "
                    f"n={evidence['getcontacts_atom_contact_count']}; "
                    f"min={evidence['getcontacts_min_distance_A']} A"
                )
            else:
                template_residue = ""
                ev_text = ""
            rows.append(
                {
                    "molecule": "C5AR",
                    "role": role,
                    "docking_chain": "X",
                    "docking_structure_resid": struct_resid,
                    "docking_structure_residue": f"X:{residue['resname']}{struct_resid}",
                    "paper_or_crystal_residue": c5ar_paper_label(str(residue["resname"]), struct_resid),
                    "paper_or_crystal_resid": struct_resid + C5AR_PAPER_NUMBER_OFFSET,
                    "source_template_residue": template_residue,
                    "evidence": ev_text,
                    "notes": "C5AR paper numbering = current C5AR-only structure resid + 33",
                }
            )

    beta_role_map = [
        ("active", HADDOCK_ACTIVE_BETA),
        ("passive", HADDOCK_PASSIVE_BETA),
    ]
    for role, beta_resids in beta_role_map:
        for beta_resid in beta_resids:
            evidence = evidence_for_beta_resid(residue_rows, beta_resid)
            if not evidence:
                continue
            residue_label = f"A:{evidence['resname']}{beta_resid}"
            ev_text = (
                f"{evidence['contact_types']}; "
                f"n={evidence['atom_contact_count_getcontacts']}; "
                f"min={evidence['min_distance_A_getcontacts']} A"
            )
            rows.append(
                {
                    "molecule": "beta-arrestin",
                    "role": role,
                    "docking_chain": "A",
                    "docking_structure_resid": beta_resid,
                    "docking_structure_residue": residue_label,
                    "paper_or_crystal_residue": residue_label,
                    "paper_or_crystal_resid": beta_resid,
                    "source_template_residue": "",
                    "evidence": ev_text,
                    "notes": "beta-arrestin uses the crystal/structure numbering in beta-GPR1 chain A",
                }
            )
    return rows


def build_report(
    *,
    residue_rows,
    pair_rows,
    mapping_rows,
    alignment_notes,
    core_c5ar_by_rep,
    all_c5ar_by_rep,
    passive_c5ar_by_rep,
    beta_core,
    beta_all,
    beta_heavy,
    haddock_numbering_rows,
) -> str:
    gpr1_all = [row for row in residue_rows if row["side"] == "gpr1"]
    gpr1_core = [row for row in gpr1_all if row["recommended_active"] == "yes"]
    beta_all_rows = [row for row in residue_rows if row["side"] == "beta_arrestin"]
    beta_core_rows = [row for row in beta_all_rows if row["recommended_active"] == "yes"]

    top_pairs = sorted(
        pair_rows,
        key=lambda row: (-int(row["atom_contact_count_getcontacts"]), float(row["min_distance_A_getcontacts"])),
    )[:12]

    lines = [
        "# β-GPR1  C5AR ",
        "",
        "## ",
        "",
        f"- `{BETA_GPR1_PDB.relative_to(ROOT)}`",
        "- chain `R` = GPR1/GPCRchain `A` = β-arrestin ",
        f"- C5AR representative`{C5AR_REPS['cluster_01_BM213_rep3'].relative_to(ROOT)}`  `{C5AR_REPS['cluster_02_C5apep_rep3'].relative_to(ROOT)}`",
        f"- getcontacts`{GETCONTACTS_TSV.relative_to(ROOT)}` `--sele 'chain R' --sele2 'chain A' --itypes all --distout`",
        "- GPR1→C5AR  chain R  C5AR chain X  BLOSUM62  C5AR docking  C5AR  β-arrestin",
        "",
        "## ",
        "",
        f"- getcontacts  β-GPR1  GPR1  `{len(gpr1_all)}`  active/core  `{len(gpr1_core)}` ",
        f"- β-arrestin  `{len(beta_all_rows)}`  active/core  `{len(beta_core_rows)}` ",
        f"- 4.5 Å heavy-atom / getcontacts  active residues",
        f"- C5AR docking PDB  chain `X`  1-281C5AR / =  + `{C5AR_PAPER_NUMBER_OFFSET}` offset  activation-core mapping  101  C5AR R134",
        "",
        "## HADDOCK  active/passive ",
        "",
        "HADDOCK  active residues  AIR HDOCK  C5AR  C5AR  docking PDB  HADDOCK ",
        "",
        "### C5AR receptor",
        "",
        "| role | C5AR  | docking PDB  |  | evidence |",
        "|---|---|---|---|---:|",
    ]

    for row in [r for r in haddock_numbering_rows if r["molecule"] == "C5AR"]:
        lines.append(
            f"| {row['role']} | {row['paper_or_crystal_residue']} | "
            f"{row['docking_structure_residue']} | {row['source_template_residue']} | "
            f"{row['evidence']} |"
        )

    lines += [
        "",
        "### β-arrestin ligand",
        "",
        "| role | β-arrestin  | HADDOCK  | evidence |",
        "|---|---|---|---:|",
    ]
    for row in [r for r in haddock_numbering_rows if r["molecule"] == "beta-arrestin"]:
        lines.append(
            f"| {row['role']} | {row['paper_or_crystal_residue']} | "
            f"{row['docking_structure_residue']} | {row['evidence']} |"
        )

    c5ar_active = [
        row for row in haddock_numbering_rows if row["molecule"] == "C5AR" and row["role"] == "active"
    ]
    c5ar_passive = [
        row for row in haddock_numbering_rows if row["molecule"] == "C5AR" and row["role"] == "passive"
    ]
    beta_active_rows = [
        row
        for row in haddock_numbering_rows
        if row["molecule"] == "beta-arrestin" and row["role"] == "active"
    ]
    beta_passive_rows = [
        row
        for row in haddock_numbering_rows
        if row["molecule"] == "beta-arrestin" and row["role"] == "passive"
    ]
    lines += [
        "",
        "###  HADDOCK residue lists",
        "",
        "- C5AR active`"
        + ", ".join(str(row["paper_or_crystal_residue"]) for row in c5ar_active)
        + "`",
        "- C5AR active`"
        + ", ".join(str(row["docking_structure_residue"]) for row in c5ar_active)
        + "`",
        "- C5AR passive`"
        + ", ".join(str(row["paper_or_crystal_residue"]) for row in c5ar_passive)
        + "`",
        "- C5AR passive`"
        + ", ".join(str(row["docking_structure_residue"]) for row in c5ar_passive)
        + "`",
        "- β-arrestin active`"
        + ", ".join(str(row["docking_structure_residue"]) for row in beta_active_rows)
        + "`",
        "- β-arrestin passive`"
        + ", ".join(str(row["docking_structure_residue"]) for row in beta_passive_rows)
        + "`",
        "",
        "## HDOCK/",
        "",
        "### C5AR receptor site residues",
        "",
    ]

    for name in C5AR_REPS:
        lines += [
            f"- `{name}` core/recommended`{format_site_list('X', core_c5ar_by_rep[name])}`",
            f"- `{name}` all getcontacts-mapped`{format_site_list('X', all_c5ar_by_rep[name])}`",
            f"- `{name}` 4.5 Å /passive`{format_site_list('X', passive_c5ar_by_rep[name])}`",
        ]

    lines += [
        "",
        "### β-arrestin ligand site residues",
        "",
        f"- core/recommended`{format_site_list('A', beta_core)}`",
        f"- all getcontacts`{format_site_list('A', beta_all)}`",
        f"- 4.5 Å /passive`{format_site_list('A', beta_heavy)}`",
        "",
        "## GPR1  C5AR ",
        "",
    ]

    for name, note in alignment_notes.items():
        lines.append(
            f"- `{name}`C5AR residues={note['residue_count']}aligned columns={note['aligned_cols']}sequence identity={note['identity']:.3f}score={note['score']:.1f}"
        )
    lines += [
        "",
        "| GPR1 getcontacts residue | evidence | C5AR/BM213 mapped | C5AR/C5apep mapped | use |",
        "|---|---:|---|---|---|",
    ]
    for row in mapping_rows:
        if row["in_getcontacts"] != "yes":
            continue
        evidence = (
            f"{row['getcontacts_contact_types']}; n={row['getcontacts_atom_contact_count']}; "
            f"min={row['getcontacts_min_distance_A']} Å"
        )
        bm = f"X:{row['cluster_01_BM213_rep3_c5ar_resname']}{row['cluster_01_BM213_rep3_c5ar_resid']}"
        c5 = f"X:{row['cluster_02_C5apep_rep3_c5ar_resname']}{row['cluster_02_C5apep_rep3_c5ar_resid']}"
        use = "active/core" if row["getcontacts_recommended_active"] == "yes" else "broader/passive"
        lines.append(f"| R:{row['gpr1_resname']}{row['gpr1_resid']} | {evidence} | {bm} | {c5} | {use} |")

    lines += [
        "",
        "## β-arrestin ",
        "",
        "| β-arrestin residue | evidence | partners | use |",
        "|---|---:|---|---|",
    ]
    for row in sorted(beta_all_rows, key=lambda item: int(item["resid"])):
        evidence = f"{row['contact_types']}; n={row['atom_contact_count_getcontacts']}; min={row['min_distance_A_getcontacts']} Å"
        use = "active/core" if row["recommended_active"] == "yes" else "broader/passive"
        lines.append(
            f"| {row['chain']}:{row['resname']}{row['resid']} | {evidence} | {row['partner_residues']} | {use} |"
        )

    lines += [
        "",
        "## ",
        "",
        "| GPR1 | β-arrestin | evidence |",
        "|---|---|---:|",
    ]
    for row in top_pairs:
        g = f"{row['gpr1_chain']}:{row['gpr1_resname']}{row['gpr1_resid']}"
        b = f"{row['beta_chain']}:{row['beta_resname']}{row['beta_resid']}"
        evidence = (
            f"{row['contact_types']}; n={row['atom_contact_count_getcontacts']}; "
            f"min={row['min_distance_A_getcontacts']} Å"
        )
        lines.append(f"| {g} | {b} | {evidence} |")

    lines += [
        "",
        "##  docking ",
        "",
        f"- β-arrestin `{(FORMAL / 'beta_arrestin_chainA_from_beta_GPR1.pdb').relative_to(ROOT)}`",
        f"- residue-level `{(OUTDIR / 'beta_gpr1_getcontacts_interface_residues.csv').relative_to(ROOT)}`",
        f"- contact-pair `{(OUTDIR / 'beta_gpr1_getcontacts_interface_pairs.csv').relative_to(ROOT)}`",
        f"- GPR1→C5AR `{(OUTDIR / 'gpr1_interface_to_c5ar_mapping.csv').relative_to(ROOT)}`",
        f"- HADDOCK active/passive `{(OUTDIR / 'haddock_active_passive_numbering.csv').relative_to(ROOT)}`",
        f"- BM213 HADDOCK receptor active `{(OUTDIR / 'haddock_active_rsite_cluster_01_BM213_rep3.txt').relative_to(ROOT)}`",
        f"- C5apep HADDOCK receptor active `{(OUTDIR / 'haddock_active_rsite_cluster_02_C5apep_rep3.txt').relative_to(ROOT)}`",
        f"- β-arrestin HADDOCK ligand active `{(OUTDIR / 'haddock_active_lsite_beta_arrestin.txt').relative_to(ROOT)}`",
        f"- BM213 HDOCK receptor core `{(OUTDIR / 'hdock_rsite_cluster_01_BM213_rep3_core.txt').relative_to(ROOT)}`",
        f"- C5apep HDOCK receptor core `{(OUTDIR / 'hdock_rsite_cluster_02_C5apep_rep3_core.txt').relative_to(ROOT)}`",
        f"- β-arrestin HDOCK ligand core `{(OUTDIR / 'hdock_lsite_beta_arrestin_core.txt').relative_to(ROOT)}`",
        "",
        "##  HDOCK/HADDOCK ",
        "",
        "- HDOCK C5AR receptor  core/recommended `rsite`β-arrestin ligand  core/recommended `lsite` top poses  all getcontacts ",
        "- HADDOCK  active residues passive residues getcontacts  AIR",
        "- BM213  C5apep  representative  C5AR  docking  β-arrestin  C5AR  intracellular pocket ",
        "-  GPR1-β-arrestin  β-arrestin-compatible docking  C5AR  β-arrestin ",
        "",
    ]
    return "\n".join(lines)


if __name__ == "__main__":
    main()
