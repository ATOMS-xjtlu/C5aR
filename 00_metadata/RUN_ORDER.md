# Recommended Run Order

Generated: 2026-07-02T20:11:02+08:00

This package is organized for manuscript submission and provenance review. It is not a one-click full rerun of all raw molecular dynamics because raw trajectories and large derived outputs are intentionally excluded.

1. Inspect MD setup provenance in `01_md_simulation_inputs/` and the system table in `00_metadata/md_system_parameters.tsv`.
2. Validate or stage the raw trajectories using `00_metadata/trajectory_manifest_3systems.tsv` and the original repository-relative paths.
3. Regenerate basic QC, RMSF, pocket, and microswitch panels with scripts in `02_basic_qc_pocket_microswitch/`.
4. Regenerate AlloViz-derived network tables and route views with scripts in `03_alloviz_networks/`.
5. Regenerate SOMMD state maps and representative selections with scripts in `04_sommd_state_analysis/`.
6. Regenerate transition-change panels for Figure 9 with scripts in `05_sommd_transition_changes/` and views in `07_pymol_views/sommd_transition_change/`.
7. Regenerate HADDOCK beta-arrestin compatibility panels with scripts in `06_haddock_beta_arrestin/`.
8. Use `07_pymol_views/` to reopen or ray-trace the structural panels.

Important: most scripts expect to be run from the repository root and use the original relative paths under `article/` and `gpcr/`.
