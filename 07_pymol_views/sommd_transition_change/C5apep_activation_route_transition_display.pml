# Display C5apep activation-route changes across formal SOM states
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
load E:/work/modeling/c5ar/md/article/analysis/sommd/pymol/formal/neuron_representatives/neuron_53_C4_C5apep_rep3_frame00018.pdb, C5apep_C4
hide everything, C5apep_C4
show cartoon, C5apep_C4
color cluster_ref, C5apep_C4
load E:/work/modeling/c5ar/md/article/analysis/sommd/pymol/formal/cluster_07_C5apep_rep2.pdb, C5apep_C7
hide everything, C5apep_C7
show cartoon, C5apep_C7
color cluster_state_a, C5apep_C7
load E:/work/modeling/c5ar/md/article/analysis/sommd/pymol/formal/cluster_03_C5apep_rep3.pdb, C5apep_C3
hide everything, C5apep_C3
show cartoon, C5apep_C3
color cluster_state_b, C5apep_C3
load E:/work/modeling/c5ar/md/article/analysis/sommd/pymol/formal/cluster_02_C5apep_rep3.pdb, C5apep_C2
hide everything, C5apep_C2
show cartoon, C5apep_C2
color cluster_state_c, C5apep_C2
load E:/work/modeling/c5ar/md/article/analysis/sommd/pymol/formal/archive_non_dominant_cluster_pdbs/cluster_05_C5apep_rep2.pdb, C5apep_C5
hide everything, C5apep_C5
show cartoon, C5apep_C5
color cluster_state_d, C5apep_C5
align C5apep_C7 and chain X and name CA, C5apep_C4 and chain X and name CA
align C5apep_C3 and chain X and name CA, C5apep_C4 and chain X and name CA
align C5apep_C2 and chain X and name CA, C5apep_C4 and chain X and name CA
align C5apep_C5 and chain X and name CA, C5apep_C4 and chain X and name CA
select C5apep_W102_all, (C5apep_C4 or C5apep_C7 or C5apep_C3 or C5apep_C2 or C5apep_C5) and chain X and resi 69 and name CA
show spheres, C5apep_W102_all
set sphere_scale, 0.74, C5apep_W102_all
color residue_hot, C5apep_W102_all
select C5apep_W102_label, C5apep_C4 and chain X and resi 69 and name CA
color residue_ref, C5apep_W102_label
label C5apep_W102_label, "W102"
select C5apep_Y222_all, (C5apep_C4 or C5apep_C7 or C5apep_C3 or C5apep_C2 or C5apep_C5) and chain X and resi 189 and name CA
show spheres, C5apep_Y222_all
set sphere_scale, 0.74, C5apep_Y222_all
color residue_hot, C5apep_Y222_all
select C5apep_Y222_label, C5apep_C4 and chain X and resi 189 and name CA
color residue_ref, C5apep_Y222_label
label C5apep_Y222_label, "Y222"
select C5apep_F251_all, (C5apep_C4 or C5apep_C7 or C5apep_C3 or C5apep_C2 or C5apep_C5) and chain X and resi 218 and name CA
show spheres, C5apep_F251_all
set sphere_scale, 0.74, C5apep_F251_all
color residue_hot, C5apep_F251_all
select C5apep_F251_label, C5apep_C4 and chain X and resi 218 and name CA
color residue_ref, C5apep_F251_label
label C5apep_F251_label, "F251"
select C5apep_R134_all, (C5apep_C4 or C5apep_C7 or C5apep_C3 or C5apep_C2 or C5apep_C5) and chain X and resi 101 and name CA
show spheres, C5apep_R134_all
set sphere_scale, 0.74, C5apep_R134_all
color residue_hot, C5apep_R134_all
select C5apep_R134_label, C5apep_C4 and chain X and resi 101 and name CA
color residue_ref, C5apep_R134_label
label C5apep_R134_label, "R134"
select C5apep_Y300_all, (C5apep_C4 or C5apep_C7 or C5apep_C3 or C5apep_C2 or C5apep_C5) and chain X and resi 267 and name CA
show spheres, C5apep_Y300_all
set sphere_scale, 0.62, C5apep_Y300_all
color residue_hot, C5apep_Y300_all
select C5apep_Y300_label, C5apep_C4 and chain X and resi 267 and name CA
color residue_ref, C5apep_Y300_label
label C5apep_Y300_label, "Y300"
select C5apep_N296_all, (C5apep_C4 or C5apep_C7 or C5apep_C3 or C5apep_C2 or C5apep_C5) and chain X and resi 263 and name CA
show spheres, C5apep_N296_all
set sphere_scale, 0.62, C5apep_N296_all
color residue_hot, C5apep_N296_all
select C5apep_N296_label, C5apep_C4 and chain X and resi 263 and name CA
color residue_ref, C5apep_N296_label
label C5apep_N296_label, "N296"
select C5apep_W255_all, (C5apep_C4 or C5apep_C7 or C5apep_C3 or C5apep_C2 or C5apep_C5) and chain X and resi 222 and name CA
show spheres, C5apep_W255_all
set sphere_scale, 0.62, C5apep_W255_all
color residue_hot, C5apep_W255_all
select C5apep_W255_label, C5apep_C4 and chain X and resi 222 and name CA
color residue_ref, C5apep_W255_label
label C5apep_W255_label, "W255"
select C5apep_Y258_all, (C5apep_C4 or C5apep_C7 or C5apep_C3 or C5apep_C2 or C5apep_C5) and chain X and resi 225 and name CA
show spheres, C5apep_Y258_all
set sphere_scale, 0.62, C5apep_Y258_all
color residue_hot, C5apep_Y258_all
select C5apep_Y258_label, C5apep_C4 and chain X and resi 225 and name CA
color residue_ref, C5apep_Y258_label
label C5apep_Y258_label, "Y258"
# W102-centered route
# W102-F251; weighted |delta|=2.718 A
distance pair_C5apep_C4_W102_F251, C5apep_C4 and chain X and resi 69 and name CA, C5apep_C4 and chain X and resi 218 and name CA
hide labels, pair_C5apep_C4_W102_F251
color route_entry, pair_C5apep_C4_W102_F251
set dash_width, 7.35, pair_C5apep_C4_W102_F251
set dash_radius, 0.147, pair_C5apep_C4_W102_F251
set dash_gap, 0.00, pair_C5apep_C4_W102_F251
# W102-W255; weighted |delta|=2.570 A
distance pair_C5apep_C4_W102_W255, C5apep_C4 and chain X and resi 69 and name CA, C5apep_C4 and chain X and resi 222 and name CA
hide labels, pair_C5apep_C4_W102_W255
color route_entry, pair_C5apep_C4_W102_W255
set dash_width, 6.63, pair_C5apep_C4_W102_W255
set dash_radius, 0.133, pair_C5apep_C4_W102_W255
set dash_gap, 0.00, pair_C5apep_C4_W102_W255
# W102-Y258; weighted |delta|=2.351 A
distance pair_C5apep_C4_W102_Y258, C5apep_C4 and chain X and resi 69 and name CA, C5apep_C4 and chain X and resi 225 and name CA
hide labels, pair_C5apep_C4_W102_Y258
color route_entry, pair_C5apep_C4_W102_Y258
set dash_width, 5.56, pair_C5apep_C4_W102_Y258
set dash_radius, 0.111, pair_C5apep_C4_W102_Y258
set dash_gap, 0.00, pair_C5apep_C4_W102_Y258
# W102-N296; weighted |delta|=2.132 A
distance pair_C5apep_C4_W102_N296, C5apep_C4 and chain X and resi 69 and name CA, C5apep_C4 and chain X and resi 263 and name CA
hide labels, pair_C5apep_C4_W102_N296
color route_entry, pair_C5apep_C4_W102_N296
set dash_width, 4.48, pair_C5apep_C4_W102_N296
set dash_radius, 0.090, pair_C5apep_C4_W102_N296
set dash_gap, 0.00, pair_C5apep_C4_W102_N296
# W102-Y300; weighted |delta|=1.880 A
distance pair_C5apep_C4_W102_Y300, C5apep_C4 and chain X and resi 69 and name CA, C5apep_C4 and chain X and resi 267 and name CA
hide labels, pair_C5apep_C4_W102_Y300
color route_entry, pair_C5apep_C4_W102_Y300
set dash_width, 3.25, pair_C5apep_C4_W102_Y300
set dash_radius, 0.065, pair_C5apep_C4_W102_Y300
set dash_gap, 0.00, pair_C5apep_C4_W102_Y300
# TM5/TM6-to-NPxxY
# R134-Y300; weighted |delta|=2.769 A
distance pair_C5apep_C4_R134_Y300, C5apep_C4 and chain X and resi 101 and name CA, C5apep_C4 and chain X and resi 267 and name CA
hide labels, pair_C5apep_C4_R134_Y300
color route_core, pair_C5apep_C4_R134_Y300
set dash_width, 7.60, pair_C5apep_C4_R134_Y300
set dash_radius, 0.152, pair_C5apep_C4_R134_Y300
set dash_gap, 0.00, pair_C5apep_C4_R134_Y300
# Y222-Y300; weighted |delta|=2.247 A
distance pair_C5apep_C4_Y222_Y300, C5apep_C4 and chain X and resi 189 and name CA, C5apep_C4 and chain X and resi 267 and name CA
hide labels, pair_C5apep_C4_Y222_Y300
color route_core, pair_C5apep_C4_Y222_Y300
set dash_width, 5.05, pair_C5apep_C4_Y222_Y300
set dash_radius, 0.101, pair_C5apep_C4_Y222_Y300
set dash_gap, 0.00, pair_C5apep_C4_Y222_Y300
# R134-N296; weighted |delta|=1.962 A
distance pair_C5apep_C4_R134_N296, C5apep_C4 and chain X and resi 101 and name CA, C5apep_C4 and chain X and resi 263 and name CA
hide labels, pair_C5apep_C4_R134_N296
color route_core, pair_C5apep_C4_R134_N296
set dash_width, 3.65, pair_C5apep_C4_R134_N296
set dash_radius, 0.073, pair_C5apep_C4_R134_N296
set dash_gap, 0.00, pair_C5apep_C4_R134_N296
# Y222-F251; weighted |delta|=1.832 A
distance pair_C5apep_C4_Y222_F251, C5apep_C4 and chain X and resi 189 and name CA, C5apep_C4 and chain X and resi 218 and name CA
hide labels, pair_C5apep_C4_Y222_F251
color route_core, pair_C5apep_C4_Y222_F251
set dash_width, 3.01, pair_C5apep_C4_Y222_F251
set dash_radius, 0.060, pair_C5apep_C4_Y222_F251
set dash_gap, 0.00, pair_C5apep_C4_Y222_F251
# TM5-to-CWxP coupling
# Y222-W255; weighted |delta|=1.706 A
distance pair_C5apep_C4_Y222_W255, C5apep_C4 and chain X and resi 189 and name CA, C5apep_C4 and chain X and resi 222 and name CA
hide labels, pair_C5apep_C4_Y222_W255
color route_bridge, pair_C5apep_C4_Y222_W255
set dash_width, 2.40, pair_C5apep_C4_Y222_W255
set dash_radius, 0.048, pair_C5apep_C4_Y222_W255
set dash_gap, 0.00, pair_C5apep_C4_Y222_W255
# Y222-Y258; weighted |delta|=1.769 A
distance pair_C5apep_C4_Y222_Y258, C5apep_C4 and chain X and resi 189 and name CA, C5apep_C4 and chain X and resi 225 and name CA
hide labels, pair_C5apep_C4_Y222_Y258
color route_bridge, pair_C5apep_C4_Y222_Y258
set dash_width, 2.71, pair_C5apep_C4_Y222_Y258
set dash_radius, 0.054, pair_C5apep_C4_Y222_Y258
set dash_gap, 0.00, pair_C5apep_C4_Y222_Y258
select activation_change_residues, C5apep_W102_all or C5apep_Y222_all or C5apep_F251_all or C5apep_R134_all or C5apep_Y300_all or C5apep_N296_all or C5apep_W255_all or C5apep_Y258_all
show spheres, activation_change_residues
zoom activation_change_residues, 8
orient activation_change_residues
# note 1: C5apep: transitions are more W102-centered.
# note 2: W102 couples simultaneously to TM6/CWxP/PIF and TM7/NPxxY residues.
# note 3: C4<->C7 is the dominant W102-linked transition.
