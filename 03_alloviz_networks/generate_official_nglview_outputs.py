#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

import AlloViz
import nbformat as nbf
import nglview as nv
import numpy as np

from alloviz_run_serial import install_serial_analysis_patch
from generate_official_views import (
    DELTA_ELEMENT,
    DELTA_FILTERING,
    DELTA_METRIC,
    DELTA_PKG,
    OUTROOT,
    ROOT,
    SYSTEMS,
    create_protein,
    safe_write,
)


OUTDIR = OUTROOT / "official_views/nglview"

COMPARISONS = [
    ("BM213", "apo"),
    ("C5apep", "apo"),
    ("BM213", "C5apep"),
]


def build_delta_view(left_key: str, right_key: str, replicate: int):
    outdir = OUTROOT / "official_views/delta" / f"{left_key}_vs_{right_key}" / f"rep{replicate}"
    left = create_protein(
        SYSTEMS[left_key],
        replicate,
        outdir / left_key,
        pkgs=[DELTA_PKG],
        filterings=[DELTA_FILTERING],
    )
    right = create_protein(
        SYSTEMS[right_key],
        replicate,
        outdir / right_key,
        pkgs=[DELTA_PKG],
        filterings=[DELTA_FILTERING],
    )
    delta = AlloViz.Delta(left, right)
    safe_write(outdir / "delta_alignment.txt", str(delta._aln) + "\n")
    view = delta.view(pkg=DELTA_PKG, metric=DELTA_METRIC, filtering=DELTA_FILTERING, element=DELTA_ELEMENT)
    return delta, view


def write_notebook(replicate: int) -> Path:
    nb = nbf.v4.new_notebook()
    nb["metadata"]["kernelspec"] = {
        "display_name": "alloviz",
        "language": "python",
        "name": "python3",
    }
    nb["metadata"]["language_info"] = {"name": "python", "pygments_lexer": "ipython3"}
    cells = [
        nbf.v4.new_markdown_cell(
            "# C5aR1 AlloViz official Delta views\n\n"
            "Run this notebook with `/mnt/e/app/conda-envs/alloviz/bin/python` as the kernel. "
            "The cells follow the official AlloViz `AlloViz.Delta(...).view(...)` style."
        ),
        nbf.v4.new_code_cell(
            "import sys\n"
            "from pathlib import Path\n"
            "sys.path.insert(0, str(Path('/mnt/e/work/modeling/c5ar/md/article/analysis/alloviz/scripts')))\n"
            "import AlloViz\n"
            "from alloviz_run_serial import install_serial_analysis_patch\n"
            "from generate_official_views import SYSTEMS, create_protein, DELTA_PKG, DELTA_FILTERING, DELTA_ELEMENT, DELTA_METRIC, OUTROOT\n"
            "install_serial_analysis_patch()\n"
        ),
    ]
    for left, right in COMPARISONS:
        comp = f"{left}_vs_{right}"
        cells.append(nbf.v4.new_markdown_cell(f"## {comp} rep{replicate}"))
        cells.append(
            nbf.v4.new_code_cell(
                f"left = create_protein(SYSTEMS['{left}'], {replicate}, OUTROOT / 'official_views/delta/{comp}/rep{replicate}/{left}', pkgs=[DELTA_PKG], filterings=[DELTA_FILTERING])\n"
                f"right = create_protein(SYSTEMS['{right}'], {replicate}, OUTROOT / 'official_views/delta/{comp}/rep{replicate}/{right}', pkgs=[DELTA_PKG], filterings=[DELTA_FILTERING])\n"
                "delta = AlloViz.Delta(left, right)\n"
                "print(delta._aln)\n"
            )
        )
        cells.append(
            nbf.v4.new_code_cell(
                "delta.view(pkg=DELTA_PKG, metric=DELTA_METRIC, filtering=DELTA_FILTERING, element=DELTA_ELEMENT)"
            )
        )
    nb["cells"] = cells
    out = OUTDIR / f"c5ar1_alloviz_delta_views_rep{replicate}.ipynb"
    OUTDIR.mkdir(parents=True, exist_ok=True)
    nbf.write(nb, out)
    return out


def main() -> int:
    install_serial_analysis_patch()
    old_default = json.JSONEncoder.default

    def json_default(self, obj):
        if isinstance(obj, np.generic):
            return obj.item()
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return old_default(self, obj)

    json.JSONEncoder.default = json_default

    OUTDIR.mkdir(parents=True, exist_ok=True)
    statuses = []
    replicate = 1
    for left, right in COMPARISONS:
        comp = f"{left}_vs_{right}"
        print(f"[nglview html] {comp} rep{replicate}", flush=True)
        html = OUTDIR / f"{comp}_rep{replicate}_delta_view.html"
        delta, view = build_delta_view(left, right, replicate)
        nv.write_html(str(html), [view])
        statuses.append(
            {
                "comparison": comp,
                "replicate": f"rep{replicate}",
                "status": "success",
                "html": str(html.relative_to(OUTROOT)),
                "alignment": str(
                    (OUTROOT / "official_views/delta" / comp / f"rep{replicate}" / "delta_alignment.txt").relative_to(OUTROOT)
                ),
                "view_repr": repr(view),
            }
        )
    notebook = write_notebook(replicate)
    manifest = {
        "note": "Standalone nglview HTML files and a tutorial-style notebook for AlloViz.Delta.view.",
        "delta_call": 'delta.view(pkg="correlationplus_CA_Pear", metric="btw", filtering="No_Sequence_Neighbors", element="edges")',
        "notebook": str(notebook.relative_to(OUTROOT)),
        "statuses": statuses,
    }
    safe_write(OUTDIR / "nglview_manifest.json", json.dumps(manifest, indent=2, ensure_ascii=False) + "\n")
    print(f"Wrote {OUTDIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
