# BM213/C5apep formal SOM activation-route transition changes

## Method

- Used the 9-cluster formal transition network edge threshold shown in the figure: directed between-cluster edges with transitions >= 5.
- Quantified state-change residues from the 91 SOM activation-core distance features by comparing system-specific cluster mean distances before and after each transition.
- Residue magnitude is the mean absolute change across all activation-core distance features involving that residue; system-level magnitude is weighted by transition count.
- PDB support was computed by C-alpha aligning the representative structures and measuring per-residue C-alpha displacement for activation-core and AlloViz-route residues.

## Dominant Transitions

### BM213
- C4 -> C7: 24 transitions, P=0.030
- C7 -> C4: 24 transitions, P=0.421
- C4 -> C5: 12 transitions, P=0.015
- C5 -> C4: 10 transitions, P=0.013
- C5 -> C1: 6 transitions, P=0.008
- C7 -> C8: 6 transitions, P=0.105
- C8 -> C7: 6 transitions, P=0.003
- C1 -> C5: 5 transitions, P=0.004

### C5apep
- C4 -> C7: 16 transitions, P=0.080
- C7 -> C4: 14 transitions, P=0.004
- C3 -> C2: 9 transitions, P=0.042
- C2 -> C3: 8 transitions, P=0.005
- C3 -> C5: 5 transitions, P=0.023

## Main Residues By System

### BM213
- W102: weighted mean |delta|=1.05 A; max transition |delta|=1.35 A; top pairs: W102-W255:-1.93A; W102-F251:-1.85A; W102-Y258:-1.74A | W102-W255:+1.93A; W102-F251:+1.85A; W102-Y258:+1.74A
- R134: weighted mean |delta|=0.81 A; max transition |delta|=0.99 A; top pairs: R134-Y300:+3.00A; R134-N296:+2.28A; R134-F251:+1.69A | R134-Y300:-3.00A; R134-N296:-2.28A; R134-F251:-1.69A
- Y300: weighted mean |delta|=0.77 A; max transition |delta|=2.24 A; top pairs: Y222-Y300:-6.02A; R134-Y300:-5.26A; F251-Y300:-3.25A | Y222-Y300:+6.02A; R134-Y300:+5.26A; F251-Y300:+3.25A
- F251: weighted mean |delta|=0.75 A; max transition |delta|=0.90 A; top pairs: F251-Y300:+1.92A; F251-N296:+1.77A; R134-F251:+1.69A | F251-Y300:-1.92A; F251-N296:-1.77A; R134-F251:-1.69A
- Y222: weighted mean |delta|=0.71 A; max transition |delta|=1.19 A; top pairs: Y222-Y300:-6.02A; Y222-N296:-2.63A; Y222-Y290:-1.11A | Y222-Y300:+6.02A; Y222-N296:+2.63A; Y222-Y290:+1.11A
- N296: weighted mean |delta|=0.64 A; max transition |delta|=1.52 A; top pairs: Y222-N296:-2.63A; W255-N296:-2.29A; F251-N296:-2.16A | Y222-N296:+2.63A; W255-N296:+2.29A; F251-N296:+2.16A
- Y258: weighted mean |delta|=0.63 A; max transition |delta|=0.80 A; top pairs: R134-Y258:+1.94A; W102-Y258:-1.74A; Y222-Y258:+1.73A | R134-Y258:-1.94A; W102-Y258:+1.74A; Y222-Y258:-1.73A
- W255: weighted mean |delta|=0.58 A; max transition |delta|=0.67 A; top pairs: W102-W255:-1.93A; R134-W255:+1.81A; Y222-W255:+1.52A | W102-W255:+1.93A; R134-W255:-1.81A; Y222-W255:-1.52A

### C5apep
- W102: weighted mean |delta|=1.89 A; max transition |delta|=2.74 A; top pairs: W102-F251:-3.84A; W102-W255:-3.51A; W102-N296:-3.10A | W102-F251:+3.84A; W102-W255:+3.51A; W102-N296:+3.10A
- Y222: weighted mean |delta|=1.21 A; max transition |delta|=1.39 A; top pairs: Y222-F251:+2.52A; Y222-Y258:+2.24A; Y222-W255:+2.20A | Y222-F251:-2.52A; Y222-Y258:-2.24A; Y222-W255:-2.20A
- F251: weighted mean |delta|=1.18 A; max transition |delta|=1.43 A; top pairs: W102-F251:-3.84A; Y222-F251:+2.52A; R134-F251:+1.46A | W102-F251:+3.84A; Y222-F251:-2.52A; R134-F251:-1.46A
- R134: weighted mean |delta|=1.10 A; max transition |delta|=1.27 A; top pairs: W102-R134:-2.27A; R134-Y300:+2.17A; R134-N296:+2.07A | W102-R134:+2.27A; R134-Y300:-2.17A; R134-N296:-2.07A
- Y300: weighted mean |delta|=1.04 A; max transition |delta|=1.34 A; top pairs: R134-Y300:-4.17A; Y222-Y300:-3.27A; F251-Y300:-2.24A | R134-Y300:+4.17A; Y222-Y300:+3.27A; F251-Y300:+2.24A
- N296: weighted mean |delta|=0.99 A; max transition |delta|=1.17 A; top pairs: R134-N296:-2.21A; Y222-N296:-1.93A; F251-N296:-1.88A | R134-N296:+2.21A; Y222-N296:+1.93A; F251-N296:+1.88A
- W255: weighted mean |delta|=0.97 A; max transition |delta|=1.12 A; top pairs: W102-W255:-3.51A; Y222-W255:+2.20A; R134-W255:+1.93A | W102-W255:+3.51A; Y222-W255:-2.20A; R134-W255:-1.93A
- Y258: weighted mean |delta|=0.86 A; max transition |delta|=1.03 A; top pairs: W102-Y258:-3.07A; Y222-Y258:+2.24A; R134-Y258:+1.97A | W102-Y258:+3.07A; Y222-Y258:-2.24A; R134-Y258:-1.97A

## Main Residue Pairs By System

### BM213
- R134-Y300: weighted |delta|=1.78 A; weighted signed delta=0.12 A (increases)
- Y222-Y300: weighted |delta|=1.54 A; weighted signed delta=0.12 A (increases)
- W102-Y300: weighted |delta|=1.37 A; weighted signed delta=-0.05 A (decreases)
- W102-N292: weighted |delta|=1.30 A; weighted signed delta=-0.03 A (decreases)
- W102-F251: weighted |delta|=1.29 A; weighted signed delta=-0.02 A (decreases)
- W102-Y290: weighted |delta|=1.27 A; weighted signed delta=-0.03 A (decreases)
- W102-N296: weighted |delta|=1.24 A; weighted signed delta=-0.02 A (decreases)
- W102-W255: weighted |delta|=1.23 A; weighted signed delta=-0.01 A (decreases)

### C5apep
- R134-Y300: weighted |delta|=2.77 A; weighted signed delta=0.01 A (increases)
- W102-F251: weighted |delta|=2.72 A; weighted signed delta=-0.12 A (decreases)
- W102-W255: weighted |delta|=2.57 A; weighted signed delta=-0.09 A (decreases)
- W102-Y258: weighted |delta|=2.35 A; weighted signed delta=-0.05 A (decreases)
- Y222-Y300: weighted |delta|=2.25 A; weighted signed delta=0.04 A (increases)
- W102-N296: weighted |delta|=2.13 A; weighted signed delta=-0.21 A (decreases)
- R134-N296: weighted |delta|=1.96 A; weighted signed delta=0.16 A (increases)
- W102-Y300: weighted |delta|=1.88 A; weighted signed delta=-0.19 A (decreases)

## Per-Transition Top Residues

### BM213
- C4 -> C7: residues W102 (1.35 A), Y258 (0.80 A), R134 (0.78 A), F251 (0.75 A), W255 (0.67 A); pairs R134-Y258 (+1.94 A), W102-W255 (-1.93 A), W102-F251 (-1.85 A), R134-W255 (+1.81 A)
- C7 -> C4: residues W102 (1.35 A), Y258 (0.80 A), R134 (0.78 A), F251 (0.75 A), W255 (0.67 A); pairs R134-Y258 (-1.94 A), W102-W255 (+1.93 A), W102-F251 (+1.85 A), R134-W255 (-1.81 A)
- C4 -> C5: residues R134 (0.99 A), Y300 (0.96 A), F251 (0.90 A), I91 (0.88 A), Y222 (0.87 A); pairs R134-Y300 (+3.00 A), Y222-Y300 (+2.72 A), R134-N296 (+2.28 A), Y222-N296 (+2.12 A)
- C5 -> C4: residues R134 (0.99 A), Y300 (0.96 A), F251 (0.90 A), I91 (0.88 A), Y222 (0.87 A); pairs R134-Y300 (-3.00 A), Y222-Y300 (-2.72 A), R134-N296 (-2.28 A), Y222-N296 (-2.12 A)
- C5 -> C1: residues Y300 (2.24 A), N296 (1.52 A), Y222 (1.19 A), W102 (0.88 A), R134 (0.76 A); pairs Y222-Y300 (+6.02 A), R134-Y300 (+5.26 A), F251-Y300 (+3.25 A), M120-Y300 (+2.85 A)
- C7 -> C8: residues I91 (0.98 A), S171 (0.83 A), W102 (0.79 A), N292 (0.74 A), R134 (0.68 A); pairs I91-Y290 (-1.89 A), I91-N292 (-1.83 A), S171-Y258 (+1.67 A), I91-S171 (+1.63 A)
- C8 -> C7: residues I91 (0.98 A), S171 (0.83 A), W102 (0.79 A), N292 (0.74 A), R134 (0.68 A); pairs I91-Y290 (+1.89 A), I91-N292 (+1.83 A), S171-Y258 (-1.67 A), I91-S171 (-1.63 A)
- C1 -> C5: residues Y300 (2.24 A), N296 (1.52 A), Y222 (1.19 A), W102 (0.88 A), R134 (0.76 A); pairs Y222-Y300 (-6.02 A), R134-Y300 (-5.26 A), F251-Y300 (-3.25 A), M120-Y300 (-2.85 A)

### C5apep
- C4 -> C7: residues W102 (2.74 A), F251 (1.43 A), Y222 (1.39 A), R134 (1.27 A), W255 (1.12 A); pairs W102-F251 (-3.84 A), W102-W255 (-3.51 A), W102-N296 (-3.10 A), W102-Y300 (-3.09 A)
- C7 -> C4: residues W102 (2.74 A), F251 (1.43 A), Y222 (1.39 A), R134 (1.27 A), W255 (1.12 A); pairs W102-F251 (+3.84 A), W102-W255 (+3.51 A), W102-N296 (+3.10 A), W102-Y300 (+3.09 A)
- C3 -> C2: residues Y300 (1.34 A), N296 (1.17 A), F251 (0.98 A), R134 (0.98 A), Y222 (0.94 A); pairs R134-Y300 (+4.17 A), Y222-Y300 (+3.27 A), F251-Y300 (+2.24 A), R134-N296 (+2.21 A)
- C2 -> C3: residues Y300 (1.34 A), N296 (1.17 A), F251 (0.98 A), R134 (0.98 A), Y222 (0.94 A); pairs R134-Y300 (-4.17 A), Y222-Y300 (-3.27 A), F251-Y300 (-2.24 A), R134-N296 (-2.21 A)
- C3 -> C5: residues Y222 (0.99 A), Y300 (0.93 A), N296 (0.77 A), S171 (0.65 A), Y290 (0.59 A); pairs R134-Y300 (-1.59 A), Y222-Y290 (+1.55 A), M120-Y300 (-1.46 A), Y222-N292 (+1.41 A)

## Representative PDB Sources

### BM213
- C1: cluster_level_system_specific - cluster_01_BM213_rep3.pdb
- C4: cluster_level_system_specific - cluster_04_BM213_rep3.pdb
- C5: cluster_level_system_specific - cluster_05_BM213_rep3.pdb
- C7: archive_non_dominant_cluster - cluster_07_BM213_rep2.pdb
- C8: cluster_level_system_specific - cluster_08_BM213_rep2.pdb

### C5apep
- C2: cluster_level_system_specific - cluster_02_C5apep_rep3.pdb
- C3: cluster_level_system_specific - cluster_03_C5apep_rep3.pdb
- C4: neuron_system_specific - neuron_53_C4_C5apep_rep3_frame00018.pdb
- C5: archive_non_dominant_cluster - cluster_05_C5apep_rep2.pdb
- C7: cluster_level_system_specific - cluster_07_C5apep_rep2.pdb

## PDB C-alpha Displacement Support

### BM213
- C1 -> C5: Y300 (4.80 A), N296 (3.58 A), V243 (2.82 A), Y222 (2.74 A), R206 (2.52 A)
- C4 -> C5: Y300 (3.90 A), N296 (2.91 A), Y258 (2.31 A), R206 (2.21 A), N292 (2.08 A)
- C4 -> C7: S171 (2.81 A), R206 (2.00 A), Y290 (1.75 A), I91 (1.60 A), W102 (1.46 A)
- C5 -> C1: Y300 (4.80 A), N296 (3.58 A), V243 (2.82 A), Y222 (2.74 A), R206 (2.52 A)
- C5 -> C4: Y300 (3.90 A), N296 (2.91 A), Y258 (2.31 A), R206 (2.21 A), N292 (2.08 A)
- C7 -> C4: S171 (2.81 A), R206 (2.00 A), Y290 (1.75 A), I91 (1.60 A), W102 (1.46 A)
- C7 -> C8: representative PDB C-alpha coordinates are identical; PDB displacement is not informative for this edge.
- C8 -> C7: representative PDB C-alpha coordinates are identical; PDB displacement is not informative for this edge.

### C5apep
- C2 -> C3: Y300 (3.61 A), I91 (2.82 A), N296 (2.69 A), R310 (2.64 A), N292 (1.80 A)
- C3 -> C2: Y300 (3.61 A), I91 (2.82 A), N296 (2.69 A), R310 (2.64 A), N292 (1.80 A)
- C3 -> C5: N296 (6.76 A), Y290 (6.70 A), R206 (6.55 A), R310 (6.51 A), I116 (6.07 A)
- C4 -> C7: R310 (2.76 A), W102 (2.40 A), I91 (2.04 A), N292 (1.91 A), R206 (1.67 A)
- C7 -> C4: R310 (2.76 A), W102 (2.40 A), I91 (2.04 A), N292 (1.91 A), R206 (1.67 A)

## Notes

- These values describe structural/feature shifts associated with SOM state transitions; they are not by themselves causal proof of allosteric signal propagation.
- The PDB displacement table is a representative-structure check. The CSV feature tables are the primary quantitative basis because they average all frames assigned to each system-specific cluster.
