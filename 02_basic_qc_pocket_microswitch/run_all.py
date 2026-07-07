from __future__ import annotations

import importlib
import sys
from pathlib import Path

from figure_utils import LOG_DIR, TABLE_DIR, expected_outputs, write_figure_manifest

PLOT_MODULES = [
    "plot_fig1_overview",
    "plot_fig2_pocket_contacts",
    "plot_fig3_microswitch",
    "plot_fig4_fes_projection",
    "plot_fig5_mdpath_network",
    "plot_figS1_qc",
    "plot_c5ar_rmsf_regions",
    "plot_figS2_hbond_timeseries",
    "plot_figS3_microswitch_replicates",
    "plot_figS4_fes_basin_definition",
    "plot_figS5_mdpath_diagnostics",
]


def main() -> None:
    log_lines = []
    import collect_data

    log_lines.append("Collecting data tables...")
    collect_data.main()

    for module_name in PLOT_MODULES:
        log_lines.append(f"Rendering {module_name}.py")
        module = importlib.import_module(module_name)
        module.main()

    write_figure_manifest()
    missing = []
    for figure, path, script, _inputs in expected_outputs():
        if not path.exists() or path.stat().st_size == 0:
            missing.append(f"{figure}: {path}")
    if missing:
        log_lines.append("FAILED outputs:")
        log_lines.extend(missing)
        (LOG_DIR / "run_all.log").write_text("\n".join(log_lines) + "\n")
        raise SystemExit(1)
    log_lines.append("All PNG outputs generated successfully.")
    (LOG_DIR / "run_all.log").write_text("\n".join(log_lines) + "\n")


if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    main()
