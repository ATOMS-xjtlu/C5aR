# C5aR 2 us Input Manifest

This folder records the input layer for later C5aR calculations requested for BM213 and c5apep.

- Manifest: `trajectory_manifest.tsv`
- Validation script: `validate_inputs.py`
- Default validation output: `input_validation.tsv`

Scope:

- Systems: BM213 and c5apep.
- Replicates: three independent dry trajectories per system.
- Trajectory length: 2000 ns each.
- BM213 trajectories use `2wframe.xtc` in each replicate directory.
- c5apep uses the 2 us `2wframe` series: `2wframe.xtc`, `2wframe2.xtc`, and `2wframe3.xtc`.
- c5apep rep3 intentionally excludes `3us.xtc` for this 2000 ns input set.

Run validation:

```bash
/mnt/e/app/conda-envs/alloviz/bin/python article/analysis/c5ar_2us_inputs/validate_inputs.py
```

The manifest is intended as the shared input table for later RMSD, RMSF, distance, hydrogen-bond, contact, or network calculations.
