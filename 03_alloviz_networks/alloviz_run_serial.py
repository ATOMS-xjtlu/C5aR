#!/usr/bin/env python3
"""Small CLI wrapper for repeatable AlloViz trajectory analysis."""

from __future__ import annotations

import argparse
import importlib
import json
import os
import sys
from datetime import datetime
from pathlib import Path
from typing import Any


def split_csv(value: str | None, default: list[str] | None = None) -> list[str]:
    if value is None or value == "":
        return default or []
    return [item.strip() for item in value.split(",") if item.strip()]


def alloviz_arg(values: list[Any]) -> str | list[Any]:
    if len(values) == 1:
        return values[0]
    return values


def parse_filterings(value: str) -> list[str | list[str]]:
    parsed: list[str | list[str]] = []
    for item in split_csv(value):
        if "+" in item:
            parsed.append([part.strip() for part in item.split("+") if part.strip()])
        else:
            parsed.append(item)
    return parsed


def filtering_attr(value: str | list[str]) -> str:
    return "_".join(value) if isinstance(value, list) else value


def load_special_res(path: str | None) -> dict[str, str] | None:
    if not path:
        return None
    with open(path, "r", encoding="utf-8") as handle:
        data = json.load(handle)
    if not isinstance(data, dict):
        raise ValueError("--special-res-json must contain a JSON object")
    return {str(k): str(v) for k, v in data.items()}


def resolve_existing(path: str) -> str:
    resolved = Path(path).expanduser().resolve()
    if not resolved.exists():
        raise FileNotFoundError(str(resolved))
    return str(resolved)


def export_table(obj: Any, outfile: Path, sort_metrics: list[str]) -> dict[str, Any]:
    import pandas as pd

    frame = pd.DataFrame(obj)
    outfile.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(outfile)

    top_file = outfile.with_name(outfile.stem + "_top20.csv")
    metric_cols = [col for col in sort_metrics if col in frame.columns]
    if metric_cols:
        frame.sort_values(metric_cols[0], ascending=False).head(20).to_csv(top_file)
    else:
        frame.head(20).to_csv(top_file)

    return {
        "file": str(outfile),
        "top20_file": str(top_file),
        "rows": int(frame.shape[0]),
        "columns": list(map(str, frame.columns)),
    }


def print_methods() -> int:
    try:
        from AlloViz.AlloViz import Analysis, utils
    except Exception as exc:  # pragma: no cover - depends on environment
        print(f"Could not import AlloViz: {exc}", file=sys.stderr)
        return 2

    print("Packages:")
    for pkg in sorted(utils.pkgsl):
        print(f"  {pkg}")
    print("\nFilterings:")
    for filtering in utils.filteringsl:
        print(f"  {filtering}")
    print("\nNode metrics:")
    for metric in Analysis.nodes_dict:
        print(f"  {metric}")
    print("\nEdge metrics:")
    for metric in Analysis.edges_dict:
        print(f"  {metric}")
    return 0


def resolve_dotted_function(path: str) -> Any:
    module_name, _, attr_name = path.rpartition(".")
    if not module_name or not attr_name:
        raise ValueError(f"Invalid function path: {path}")
    module = importlib.import_module(module_name)
    return getattr(module, attr_name)


def install_serial_analysis_patch() -> None:
    """Avoid AlloViz's inner ProcessPoolExecutor fork for memory-limited WSL runs."""
    import pandas as pd
    from functools import partial
    from AlloViz.AlloViz import Analysis

    def serial_single_analysis(graphs, metricd, metric, elem, pq):
        metricf = partial(resolve_dotted_function(metricd["function"]), **metricd["arguments"])
        get_colname = lambda metric_name, col: f"{metric_name}_{col}" if len(graphs) > 1 else metric_name

        results = []
        failures = []
        for col, graph in graphs.items():
            colname = get_colname(metric, col)
            result = Analysis.analyze_graph((graph, metricf, colname))
            if isinstance(result, pd.Series):
                results.append(result)
            else:
                failures.append((colname, result))

        for colname, error in failures:
            print("ERROR:", pq, colname, elem, metric, "\n", error, file=sys.stderr)

        out = pd.concat(results, axis=1) if results else pd.DataFrame()
        out.to_parquet(pq)

    Analysis.single_analysis = serial_single_analysis


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run a focused AlloViz network analysis and export CSV tables."
    )
    parser.add_argument("--list-methods", action="store_true", help="Print available packages, filters, and metrics.")
    parser.add_argument("--pdb", help="Input PDB/topology file.")
    parser.add_argument("--trajs", nargs="+", help="One or more trajectory files.")
    parser.add_argument("--out", default="alloviz_out", help="Output directory.")
    parser.add_argument("--name", default=None, help="System name used by AlloViz.")
    parser.add_argument("--pkgs", default="GetContacts", help="Comma-separated AlloViz package names, or all.")
    parser.add_argument("--filterings", default="Spatially_distant", help="Comma-separated filters. Use A+B for filter combinations.")
    parser.add_argument("--elements", default="nodes,edges", help="nodes, edges, or nodes,edges.")
    parser.add_argument("--metrics", default="btw", help="btw, cfb, or all.")
    parser.add_argument("--cores", type=int, default=1, help="AlloViz parallel cores.")
    parser.add_argument("--taskcpus", type=int, default=None, help="Cores per package task for supported methods.")
    parser.add_argument("--stride", type=int, default=None, help="Trajectory stride for supported expensive methods.")
    parser.add_argument("--protein-sel", default=None, help="MDAnalysis selection for protein atoms.")
    parser.add_argument("--psf", default=None, help="Optional PSF for methods needing force-field topology.")
    parser.add_argument("--parameters", default=None, help="Optional NAMD parameter file for methods needing force-field parameters.")
    parser.add_argument("--special-res-json", default=None, help="JSON map of nonstandard residue names to one-letter codes.")
    parser.add_argument("--getcontacts-threshold", type=float, default=0.0)
    parser.add_argument("--sequence-neighbor-distance", type=int, default=5)
    parser.add_argument("--interresidue-distance", type=float, default=10.0)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if args.list_methods:
        return print_methods()

    if not args.pdb or not args.trajs:
        raise SystemExit("--pdb and --trajs are required unless --list-methods is used")

    try:
        import AlloViz
        from AlloViz.AlloViz import utils
    except Exception as exc:
        print(f"Could not import AlloViz in this Python environment: {exc}", file=sys.stderr)
        print("Use /mnt/e/app/conda-envs/alloviz/bin/python to run this script.", file=sys.stderr)
        return 2

    install_serial_analysis_patch()

    pdb = resolve_existing(args.pdb)
    trajs = [resolve_existing(path) for path in args.trajs]
    outdir = Path(args.out).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    pkg_values = split_csv(args.pkgs)
    pkgs_for_alloviz: str | list[str]
    export_pkgs: list[str]
    if len(pkg_values) == 1 and pkg_values[0].lower() == "all":
        pkgs_for_alloviz = "all"
        export_pkgs = sorted(utils.pkgsl)
    else:
        export_pkgs = [utils.pkgname(pkg) for pkg in pkg_values]
        pkgs_for_alloviz = alloviz_arg(export_pkgs)

    filter_values = parse_filterings(args.filterings)
    filter_arg = alloviz_arg(filter_values)
    filter_attrs = [filtering_attr(item) for item in filter_values]

    element_values = split_csv(args.elements)
    if len(element_values) == 1 and element_values[0].lower() == "all":
        element_values = ["nodes", "edges"]
    element_values = [elem.lower() for elem in element_values]

    metric_values = split_csv(args.metrics)
    metrics_for_alloviz: str | list[str] = "all" if metric_values == ["all"] else alloviz_arg(metric_values)
    sort_metrics = ["btw", "cfb"] if metric_values == ["all"] else metric_values

    protein_kwargs: dict[str, Any] = {}
    if args.protein_sel:
        protein_kwargs["protein_sel"] = args.protein_sel
    if args.psf or args.parameters:
        if not (args.psf and args.parameters):
            raise SystemExit("--psf and --parameters must be provided together")
        protein_kwargs["psf"] = resolve_existing(args.psf)
        protein_kwargs["parameters"] = resolve_existing(args.parameters)
    special_res = load_special_res(args.special_res_json)
    if special_res:
        protein_kwargs["special_res"] = special_res

    calculate_kwargs: dict[str, Any] = {}
    if args.taskcpus is not None:
        calculate_kwargs["taskcpus"] = args.taskcpus
    if args.stride is not None:
        calculate_kwargs["stride"] = args.stride

    system = AlloViz.Protein(
        pdb=pdb,
        trajs=trajs,
        name=args.name,
        path=str(outdir),
        **protein_kwargs,
    )

    system.calculate(pkgs=pkgs_for_alloviz, cores=args.cores, **calculate_kwargs)
    system.filter(
        pkgs=pkgs_for_alloviz,
        filterings=filter_arg,
        GetContacts_threshold=args.getcontacts_threshold,
        Sequence_Neighbor_distance=args.sequence_neighbor_distance,
        Interresidue_distance=args.interresidue_distance,
    )
    system.analyze(
        pkgs=pkgs_for_alloviz,
        filterings=filter_arg,
        elements=alloviz_arg(element_values),
        metrics=metrics_for_alloviz,
        cores=args.cores,
    )

    exports: list[dict[str, Any]] = []
    export_dir = outdir / "exports"
    for pkg in export_pkgs:
        pkg_obj = getattr(system, pkg, None)
        if pkg_obj is None:
            continue
        for filtering in filter_attrs:
            filtering_obj = getattr(pkg_obj, filtering, None)
            if filtering_obj is None:
                continue
            for element in element_values:
                obj = getattr(filtering_obj, element, None)
                if obj is None:
                    continue
                outfile = export_dir / f"{pkg}__{filtering}__{element}.csv"
                meta = export_table(obj, outfile, sort_metrics)
                meta.update({"pkg": pkg, "filtering": filtering, "element": element})
                exports.append(meta)

    summary = {
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "pdb": pdb,
        "trajs": trajs,
        "outdir": str(outdir),
        "pkgs": export_pkgs,
        "filterings": filter_attrs,
        "elements": element_values,
        "metrics": metric_values,
        "exports": exports,
    }
    summary_path = outdir / "run_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
