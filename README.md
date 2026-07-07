# C5aR1 JCIM Figure Code Submission Package

Generated: 2026-07-02T20:11:02+08:00

This directory is the submission-oriented version of `result/figure/code`. It groups figure-generation scripts, PyMOL view files, MD input/parameter provenance, and metadata tables for the current C5aR1 manuscript figures.

## What Is Included

- `00_metadata/`: source-to-figure mapping, three-system trajectory manifest, MD system parameter table, script inventory, and recommended run order.
- `01_md_simulation_inputs/`: Amber `.mdin` files, run scripts, cpptraj/DSSP/distance/BSA `.in` files, CHARMM-GUI setup provenance, topology/coordinate files, and restraint files for apo, BM213, and C5apep.
- `02_basic_qc_pocket_microswitch/`: basic QC, RMSF, pocket, and microswitch plotting scripts.
- `03_alloviz_networks/`: AlloViz post-processing, NetworkX route analysis, and visualization-generation scripts.
- `04_sommd_state_analysis/`: SOMMD feature preparation, SOM training, summarization, population plots, and representative-selection scripts.
- `05_sommd_transition_changes/`: transition-change analysis scripts for Figure 9.
- `06_haddock_beta_arrestin/`: HADDOCK scoring, beta-arrestin LRMSD, and plotting scripts.
- `07_pymol_views/`: PyMOL view files used for structural panels.

## What Is Not Included

Raw trajectories (`.xtc`, `.trj`, `.nc`, `.dcd`) and large derived outputs are not copied here. They should be archived separately or restored at the paths recorded in `00_metadata/trajectory_manifest_3systems.tsv`.

## Submission Caveat

The code package is now suitable as a provenance-oriented submission supplement, but a fully external one-command rerun still requires the raw trajectories, the same relative directory layout, and the software environments listed in `00_metadata/RUN_ORDER.md`.
