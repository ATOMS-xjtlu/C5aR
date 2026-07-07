# Display BM213 activation-route changes across formal SOM states
# Lines connect activation-core residue pairs with largest transition-associated distance changes.
# Line thickness scales with transition-weighted mean absolute distance change.
reinitialize
bg_color white
set ray_opaque_background, off
set cartoon_fancy_helices, 1
set cartoon_transparency, 0.58
set sphere_quality, 2
set label_size, 20
set dash_round_ends, 1
set_color cluster_ref, [0.28, 0.28, 0.28]
set_color cluster_state_a, [0.55, 0.65, 0.75]
set_color cluster_state_b, [0.72, 0.72, 0.72]
set_color cluster_state_c, [0.64, 0.78, 0.86]
set_color cluster_state_d, [0.82, 0.82, 0.82]
set_color route_entry, [0.00, 0.62, 0.72]
set_color route_core, [0.80, 0.18, 0.15]
set_color route_bridge, [0.92, 0.55, 0.12]
set_color residue_hot, [1.00, 0.88, 0.18]
set_color residue_ref, [0.10, 0.10, 0.10]
load E:/work/modeling/c5ar/md/article/analysis/sommd/pymol/formal/cluster_04_BM213_rep3.pdb, BM213_C4
hide everything, BM213_C4
show cartoon, BM213_C4
color cluster_ref, BM213_C4
load E:/work/modeling/c5ar/md/article/analysis/sommd/pymol/formal/archive_non_dominant_cluster_pdbs/cluster_07_BM213_rep2.pdb, BM213_C7
hide everything, BM213_C7
show cartoon, BM213_C7
color cluster_state_a, BM213_C7
load E:/work/modeling/c5ar/md/article/analysis/sommd/pymol/formal/cluster_05_BM213_rep3.pdb, BM213_C5
hide everything, BM213_C5
show cartoon, BM213_C5
color cluster_state_b, BM213_C5
load E:/work/modeling/c5ar/md/article/analysis/sommd/pymol/formal/cluster_01_BM213_rep3.pdb, BM213_C1
hide everything, BM213_C1
show cartoon, BM213_C1
color cluster_state_c, BM213_C1
load E:/work/modeling/c5ar/md/article/analysis/sommd/pymol/formal/cluster_08_BM213_rep2.pdb, BM213_C8
hide everything, BM213_C8
show cartoon, BM213_C8
color cluster_state_d, BM213_C8
align BM213_C7 and chain X and name CA, BM213_C4 and chain X and name CA
align BM213_C5 and chain X and name CA, BM213_C4 and chain X and name CA
align BM213_C1 and chain X and name CA, BM213_C4 and chain X and name CA
align BM213_C8 and chain X and name CA, BM213_C4 and chain X and name CA
select BM213_W102_all, (BM213_C4 or BM213_C7 or BM213_C5 or BM213_C1 or BM213_C8) and chain X and resi 69 and name CA
show spheres, BM213_W102_all
set sphere_scale, 0.74, BM213_W102_all
color residue_hot, BM213_W102_all
select BM213_W102_label, BM213_C4 and chain X and resi 69 and name CA
color residue_ref, BM213_W102_label
label BM213_W102_label, "W102"
select BM213_R134_all, (BM213_C4 or BM213_C7 or BM213_C5 or BM213_C1 or BM213_C8) and chain X and resi 101 and name CA
show spheres, BM213_R134_all
set sphere_scale, 0.74, BM213_R134_all
color residue_hot, BM213_R134_all
select BM213_R134_label, BM213_C4 and chain X and resi 101 and name CA
color residue_ref, BM213_R134_label
label BM213_R134_label, "R134"
select BM213_Y300_all, (BM213_C4 or BM213_C7 or BM213_C5 or BM213_C1 or BM213_C8) and chain X and resi 267 and name CA
show spheres, BM213_Y300_all
set sphere_scale, 0.74, BM213_Y300_all
color residue_hot, BM213_Y300_all
select BM213_Y300_label, BM213_C4 and chain X and resi 267 and name CA
color residue_ref, BM213_Y300_label
label BM213_Y300_label, "Y300"
select BM213_F251_all, (BM213_C4 or BM213_C7 or BM213_C5 or BM213_C1 or BM213_C8) and chain X and resi 218 and name CA
show spheres, BM213_F251_all
set sphere_scale, 0.74, BM213_F251_all
color residue_hot, BM213_F251_all
select BM213_F251_label, BM213_C4 and chain X and resi 218 and name CA
color residue_ref, BM213_F251_label
label BM213_F251_label, "F251"
select BM213_Y222_all, (BM213_C4 or BM213_C7 or BM213_C5 or BM213_C1 or BM213_C8) and chain X and resi 189 and name CA
show spheres, BM213_Y222_all
set sphere_scale, 0.62, BM213_Y222_all
color residue_hot, BM213_Y222_all
select BM213_Y222_label, BM213_C4 and chain X and resi 189 and name CA
color residue_ref, BM213_Y222_label
label BM213_Y222_label, "Y222"
select BM213_N296_all, (BM213_C4 or BM213_C7 or BM213_C5 or BM213_C1 or BM213_C8) and chain X and resi 263 and name CA
show spheres, BM213_N296_all
set sphere_scale, 0.62, BM213_N296_all
color residue_hot, BM213_N296_all
select BM213_N296_label, BM213_C4 and chain X and resi 263 and name CA
color residue_ref, BM213_N296_label
label BM213_N296_label, "N296"
select BM213_Y258_all, (BM213_C4 or BM213_C7 or BM213_C5 or BM213_C1 or BM213_C8) and chain X and resi 225 and name CA
show spheres, BM213_Y258_all
set sphere_scale, 0.62, BM213_Y258_all
color residue_hot, BM213_Y258_all
select BM213_Y258_label, BM213_C4 and chain X and resi 225 and name CA
color residue_ref, BM213_Y258_label
label BM213_Y258_label, "Y258"
select BM213_W255_all, (BM213_C4 or BM213_C7 or BM213_C5 or BM213_C1 or BM213_C8) and chain X and resi 222 and name CA
show spheres, BM213_W255_all
set sphere_scale, 0.62, BM213_W255_all
color residue_hot, BM213_W255_all
select BM213_W255_label, BM213_C4 and chain X and resi 222 and name CA
color residue_ref, BM213_W255_label
label BM213_W255_label, "W255"
# Pocket-to-TM6/CWxP
# W102-W255; weighted |delta|=1.231 A
distance pair_BM213_C4_W102_W255, BM213_C4 and chain X and resi 69 and name CA, BM213_C4 and chain X and resi 222 and name CA
hide labels, pair_BM213_C4_W102_W255
color route_entry, pair_BM213_C4_W102_W255
set dash_width, 3.95, pair_BM213_C4_W102_W255
set dash_radius, 0.079, pair_BM213_C4_W102_W255
set dash_gap, 0.00, pair_BM213_C4_W102_W255
# W102-F251; weighted |delta|=1.293 A
distance pair_BM213_C4_W102_F251, BM213_C4 and chain X and resi 69 and name CA, BM213_C4 and chain X and resi 218 and name CA
hide labels, pair_BM213_C4_W102_F251
color route_entry, pair_BM213_C4_W102_F251
set dash_width, 4.36, pair_BM213_C4_W102_F251
set dash_radius, 0.087, pair_BM213_C4_W102_F251
set dash_gap, 0.00, pair_BM213_C4_W102_F251
# W102-Y258; weighted |delta|=1.098 A
distance pair_BM213_C4_W102_Y258, BM213_C4 and chain X and resi 69 and name CA, BM213_C4 and chain X and resi 225 and name CA
hide labels, pair_BM213_C4_W102_Y258
color route_entry, pair_BM213_C4_W102_Y258
set dash_width, 3.05, pair_BM213_C4_W102_Y258
set dash_radius, 0.061, pair_BM213_C4_W102_Y258
set dash_gap, 0.00, pair_BM213_C4_W102_Y258
# DRY/TM5-to-NPxxY
# R134-Y300; weighted |delta|=1.775 A
distance pair_BM213_C4_R134_Y300, BM213_C4 and chain X and resi 101 and name CA, BM213_C4 and chain X and resi 267 and name CA
hide labels, pair_BM213_C4_R134_Y300
color route_core, pair_BM213_C4_R134_Y300
set dash_width, 7.60, pair_BM213_C4_R134_Y300
set dash_radius, 0.152, pair_BM213_C4_R134_Y300
set dash_gap, 0.00, pair_BM213_C4_R134_Y300
# Y222-Y300; weighted |delta|=1.543 A
distance pair_BM213_C4_Y222_Y300, BM213_C4 and chain X and resi 189 and name CA, BM213_C4 and chain X and resi 267 and name CA
hide labels, pair_BM213_C4_Y222_Y300
color route_core, pair_BM213_C4_Y222_Y300
set dash_width, 6.04, pair_BM213_C4_Y222_Y300
set dash_radius, 0.121, pair_BM213_C4_Y222_Y300
set dash_gap, 0.00, pair_BM213_C4_Y222_Y300
# R134-N296; weighted |delta|=1.110 A
distance pair_BM213_C4_R134_N296, BM213_C4 and chain X and resi 101 and name CA, BM213_C4 and chain X and resi 263 and name CA
hide labels, pair_BM213_C4_R134_N296
color route_core, pair_BM213_C4_R134_N296
set dash_width, 3.13, pair_BM213_C4_R134_N296
set dash_radius, 0.063, pair_BM213_C4_R134_N296
set dash_gap, 0.00, pair_BM213_C4_R134_N296
# Y222-N296; weighted |delta|=1.001 A
distance pair_BM213_C4_Y222_N296, BM213_C4 and chain X and resi 189 and name CA, BM213_C4 and chain X and resi 263 and name CA
hide labels, pair_BM213_C4_Y222_N296
color route_core, pair_BM213_C4_Y222_N296
set dash_width, 2.40, pair_BM213_C4_Y222_N296
set dash_radius, 0.048, pair_BM213_C4_Y222_N296
set dash_gap, 0.00, pair_BM213_C4_Y222_N296
# DRY-to-TM6 coupling
# R134-Y258; weighted |delta|=1.185 A
distance pair_BM213_C4_R134_Y258, BM213_C4 and chain X and resi 101 and name CA, BM213_C4 and chain X and resi 225 and name CA
hide labels, pair_BM213_C4_R134_Y258
color route_bridge, pair_BM213_C4_R134_Y258
set dash_width, 3.64, pair_BM213_C4_R134_Y258
set dash_radius, 0.073, pair_BM213_C4_R134_Y258
set dash_gap, 0.00, pair_BM213_C4_R134_Y258
# R134-W255; weighted |delta|=1.172 A
distance pair_BM213_C4_R134_W255, BM213_C4 and chain X and resi 101 and name CA, BM213_C4 and chain X and resi 222 and name CA
hide labels, pair_BM213_C4_R134_W255
color route_bridge, pair_BM213_C4_R134_W255
set dash_width, 3.55, pair_BM213_C4_R134_W255
set dash_radius, 0.071, pair_BM213_C4_R134_W255
set dash_gap, 0.00, pair_BM213_C4_R134_W255
select activation_change_residues, BM213_W102_all or BM213_R134_all or BM213_Y300_all or BM213_F251_all or BM213_Y222_all or BM213_N296_all or BM213_Y258_all or BM213_W255_all
show spheres, activation_change_residues
zoom activation_change_residues, 8
orient activation_change_residues
# note 1: BM213: transitions emphasize a segmented route.
# note 2: Pocket/ECL1 W102 couples to TM6/CWxP/PIF, then DRY-like/TM5 couples to NPxxY.
# note 3: C5<->C1 is the strongest NPxxY-end switch.
