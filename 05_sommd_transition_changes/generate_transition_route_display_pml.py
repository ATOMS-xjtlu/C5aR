#!/usr/bin/env python3
"""Generate PyMOL display scripts for BM213/C5apep activation-route changes."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


BASE = Path(__file__).resolve().parent
ROOT = Path("/mnt/e/work/modeling/c5ar/md")

RESID = {
    "I91": 58,
    "W102": 69,
    "I116": 83,
    "M120": 87,
    "R134": 101,
    "S171": 138,
    "Y222": 189,
    "F251": 218,
    "W255": 222,
    "Y258": 225,
    "Y290": 257,
    "N292": 259,
    "N296": 263,
    "Y300": 267,
}

SYSTEMS = {
    "BM213": {
        "reference": 4,
        "clusters": [4, 7, 5, 1, 8],
        "cartoon_colors": {
            4: "cluster_ref",
            7: "cluster_state_a",
            5: "cluster_state_b",
            1: "cluster_state_c",
            8: "cluster_state_d",
        },
        "residues": ["W102", "R134", "Y300", "F251", "Y222", "N296", "Y258", "W255"],
        "pair_groups": [
            ("Pocket-to-TM6/CWxP", ["W102-W255", "W102-F251", "W102-Y258"], "route_entry"),
            ("DRY/TM5-to-NPxxY", ["R134-Y300", "Y222-Y300", "R134-N296", "Y222-N296"], "route_core"),
            ("DRY-to-TM6 coupling", ["R134-Y258", "R134-W255"], "route_bridge"),
        ],
        "notes": [
            "BM213: transitions emphasize a segmented route.",
            "Pocket/ECL1 W102 couples to TM6/CWxP/PIF, then DRY-like/TM5 couples to NPxxY.",
            "C5<->C1 is the strongest NPxxY-end switch.",
        ],
    },
    "C5apep": {
        "reference": 4,
        "clusters": [4, 7, 3, 2, 5],
        "cartoon_colors": {
            4: "cluster_ref",
            7: "cluster_state_a",
            3: "cluster_state_b",
            2: "cluster_state_c",
            5: "cluster_state_d",
        },
        "residues": ["W102", "Y222", "F251", "R134", "Y300", "N296", "W255", "Y258"],
        "pair_groups": [
            ("W102-centered route", ["W102-F251", "W102-W255", "W102-Y258", "W102-N296", "W102-Y300"], "route_entry"),
            ("TM5/TM6-to-NPxxY", ["R134-Y300", "Y222-Y300", "R134-N296", "Y222-F251"], "route_core"),
            ("TM5-to-CWxP coupling", ["Y222-W255", "Y222-Y258"], "route_bridge"),
        ],
        "notes": [
            "C5apep: transitions are more W102-centered.",
            "W102 couples simultaneously to TM6/CWxP/PIF and TM7/NPxxY residues.",
            "C4<->C7 is the dominant W102-linked transition.",
        ],
    },
}


def pml_quote(text: str) -> str:
    return text.replace("\\", "\\\\").replace('"', '\\"')


def width_from_delta(delta: float, min_delta: float, max_delta: float) -> float:
    if max_delta <= min_delta:
        return 4.0
    return 2.4 + (delta - min_delta) / (max_delta - min_delta) * 5.2


def load_inventory() -> pd.DataFrame:
    inv = pd.read_csv(BASE / "pdb_representative_inventory.csv")
    inv["pdb_path"] = inv["pdb_path"].astype(str)
    return inv


def load_pair_scores() -> pd.DataFrame:
    pairs = pd.read_csv(BASE / "system_activation_pair_change_summary.csv")
    return pairs.set_index(["system", "residue_pair"])


def cluster_path(inv: pd.DataFrame, system: str, cluster: int) -> str:
    row = inv[(inv["system"] == system) & (inv["cluster"] == cluster)]
    if row.empty:
        raise ValueError(f"Missing PDB representative for {system} C{cluster}")
    return row.iloc[0]["pdb_path"]


def pymol_path(path: str) -> str:
    """Use Windows-friendly paths because PyMOL is commonly launched on Windows."""

    if path.startswith("/mnt/e/"):
        return "E:/" + path[len("/mnt/e/") :]
    return path


def add_distance_line(lines: list[str], obj: str, pair: str, color: str, width: float, score: float) -> None:
    residue_a, residue_b = pair.split("-")
    resid_a = RESID[residue_a]
    resid_b = RESID[residue_b]
    name = f"pair_{obj}_{residue_a}_{residue_b}"
    lines.extend(
        [
            f"# {pair}; weighted |delta|={score:.3f} A",
            f"distance {name}, {obj} and chain X and resi {resid_a} and name CA, {obj} and chain X and resi {resid_b} and name CA",
            f"hide labels, {name}",
            f"color {color}, {name}",
            f"set dash_width, {width:.2f}, {name}",
            f"set dash_radius, {width / 50:.3f}, {name}",
            f"set dash_gap, 0.00, {name}",
        ]
    )


def write_system_pml(system: str, config: dict, inv: pd.DataFrame, pair_scores: pd.DataFrame) -> None:
    reference_cluster = config["reference"]
    reference_obj = f"{system}_C{reference_cluster}"
    clusters = config["clusters"]
    objects = [f"{system}_C{cluster}" for cluster in clusters]

    selected_pairs = []
    for _, pairs, _ in config["pair_groups"]:
        selected_pairs.extend(pairs)
    deltas = [
        float(pair_scores.loc[(system, pair), "weighted_abs_delta_angstrom"])
        for pair in selected_pairs
        if (system, pair) in pair_scores.index
    ]
    min_delta, max_delta = min(deltas), max(deltas)

    lines = [
        f"# Display {system} activation-route changes across formal SOM states",
        "# Lines connect activation-core residue pairs with largest transition-associated distance changes.",
        "# Line thickness scales with transition-weighted mean absolute distance change.",
        "reinitialize",
        "bg_color white",
        "set ray_opaque_background, off",
        "set cartoon_fancy_helices, 1",
        "set cartoon_transparency, 0.58",
        "set sphere_quality, 2",
        "set label_size, 20",
        "set dash_round_ends, 1",
        "set_color cluster_ref, [0.28, 0.28, 0.28]",
        "set_color cluster_state_a, [0.55, 0.65, 0.75]",
        "set_color cluster_state_b, [0.72, 0.72, 0.72]",
        "set_color cluster_state_c, [0.64, 0.78, 0.86]",
        "set_color cluster_state_d, [0.82, 0.82, 0.82]",
        "set_color route_entry, [0.00, 0.62, 0.72]",
        "set_color route_core, [0.80, 0.18, 0.15]",
        "set_color route_bridge, [0.92, 0.55, 0.12]",
        "set_color residue_hot, [1.00, 0.88, 0.18]",
        "set_color residue_ref, [0.10, 0.10, 0.10]",
    ]
    for cluster in clusters:
        obj = f"{system}_C{cluster}"
        color = config["cartoon_colors"][cluster]
        path = cluster_path(inv, system, cluster)
        lines.extend(
            [
                f"load {pymol_path(path)}, {obj}",
                f"hide everything, {obj}",
                f"show cartoon, {obj}",
                f"color {color}, {obj}",
            ]
        )
    for obj in objects:
        if obj != reference_obj:
            lines.append(f"align {obj} and chain X and name CA, {reference_obj} and chain X and name CA")

    all_obj_selection = " or ".join(objects)
    for rank, residue in enumerate(config["residues"], start=1):
        resid = RESID[residue]
        scale = 0.74 if rank <= 4 else 0.62
        selection = f"{system}_{residue}_all"
        lines.extend(
            [
                f"select {selection}, ({all_obj_selection}) and chain X and resi {resid} and name CA",
                f"show spheres, {selection}",
                f"set sphere_scale, {scale:.2f}, {selection}",
                f"color residue_hot, {selection}",
                f"select {system}_{residue}_label, {reference_obj} and chain X and resi {resid} and name CA",
                f"color residue_ref, {system}_{residue}_label",
                f"label {system}_{residue}_label, \"{residue}\"",
            ]
        )

    for group_label, pairs, color in config["pair_groups"]:
        lines.append(f"# {group_label}")
        for pair in pairs:
            if (system, pair) not in pair_scores.index:
                continue
            score = float(pair_scores.loc[(system, pair), "weighted_abs_delta_angstrom"])
            add_distance_line(
                lines,
                reference_obj,
                pair,
                color,
                width_from_delta(score, min_delta, max_delta),
                score,
            )

    lines.extend(
        [
            "select activation_change_residues, " + " or ".join(f"{system}_{residue}_all" for residue in config["residues"]),
            "show spheres, activation_change_residues",
            "zoom activation_change_residues, 8",
            "orient activation_change_residues",
        ]
    )
    for idx, note in enumerate(config["notes"], start=1):
        lines.append(f"# note {idx}: {pml_quote(note)}")

    out = BASE / f"{system}_activation_route_transition_display.pml"
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(out)


def main() -> None:
    inv = load_inventory()
    pair_scores = load_pair_scores()
    for system, config in SYSTEMS.items():
        write_system_pml(system, config, inv, pair_scores)


if __name__ == "__main__":
    main()
