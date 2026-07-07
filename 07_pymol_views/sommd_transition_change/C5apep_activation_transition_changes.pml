# C5apep major activation-core transition changes
# Spheres highlight residues with largest transition-weighted changes.
reinitialize
bg_color white
set cartoon_transparency, 0.25
set sphere_quality, 2
set label_size, 18
load /mnt/e/work/modeling/c5ar/md/article/analysis/sommd/pymol/formal/cluster_02_C5apep_rep3.pdb, C5apep_C2
hide everything, C5apep_C2
show cartoon, C5apep_C2
color orange, C5apep_C2
load /mnt/e/work/modeling/c5ar/md/article/analysis/sommd/pymol/formal/cluster_03_C5apep_rep3.pdb, C5apep_C3
hide everything, C5apep_C3
show cartoon, C5apep_C3
color orange, C5apep_C3
load /mnt/e/work/modeling/c5ar/md/article/analysis/sommd/pymol/formal/neuron_representatives/neuron_53_C4_C5apep_rep3_frame00018.pdb, C5apep_C4
hide everything, C5apep_C4
show cartoon, C5apep_C4
color orange, C5apep_C4
load /mnt/e/work/modeling/c5ar/md/article/analysis/sommd/pymol/formal/archive_non_dominant_cluster_pdbs/cluster_05_C5apep_rep2.pdb, C5apep_C5
hide everything, C5apep_C5
show cartoon, C5apep_C5
color orange, C5apep_C5
load /mnt/e/work/modeling/c5ar/md/article/analysis/sommd/pymol/formal/cluster_07_C5apep_rep2.pdb, C5apep_C7
hide everything, C5apep_C7
show cartoon, C5apep_C7
color orange, C5apep_C7
align C5apep_C3 and chain X and name CA, C5apep_C2 and chain X and name CA
align C5apep_C4 and chain X and name CA, C5apep_C2 and chain X and name CA
align C5apep_C5 and chain X and name CA, C5apep_C2 and chain X and name CA
align C5apep_C7 and chain X and name CA, C5apep_C2 and chain X and name CA
select C5apep_W102, C5apep_C2 or C5apep_C3 or C5apep_C4 or C5apep_C5 or C5apep_C7 and chain X and resi 69 and name CA
show spheres, C5apep_W102
set sphere_scale, 0.80, C5apep_W102
color red, C5apep_W102
label C5apep_W102, "W102"
select C5apep_Y222, C5apep_C2 or C5apep_C3 or C5apep_C4 or C5apep_C5 or C5apep_C7 and chain X and resi 189 and name CA
show spheres, C5apep_Y222
set sphere_scale, 0.75, C5apep_Y222
color orange, C5apep_Y222
label C5apep_Y222, "Y222"
select C5apep_F251, C5apep_C2 or C5apep_C3 or C5apep_C4 or C5apep_C5 or C5apep_C7 and chain X and resi 218 and name CA
show spheres, C5apep_F251
set sphere_scale, 0.74, C5apep_F251
color yelloworange, C5apep_F251
label C5apep_F251, "F251"
select C5apep_R134, C5apep_C2 or C5apep_C3 or C5apep_C4 or C5apep_C5 or C5apep_C7 and chain X and resi 101 and name CA
show spheres, C5apep_R134
set sphere_scale, 0.72, C5apep_R134
color tv_yellow, C5apep_R134
label C5apep_R134, "R134"
select C5apep_Y300, C5apep_C2 or C5apep_C3 or C5apep_C4 or C5apep_C5 or C5apep_C7 and chain X and resi 267 and name CA
show spheres, C5apep_Y300
set sphere_scale, 0.70, C5apep_Y300
color salmon, C5apep_Y300
label C5apep_Y300, "Y300"
select C5apep_N296, C5apep_C2 or C5apep_C3 or C5apep_C4 or C5apep_C5 or C5apep_C7 and chain X and resi 263 and name CA
show spheres, C5apep_N296
set sphere_scale, 0.69, C5apep_N296
color raspberry, C5apep_N296
label C5apep_N296, "N296"
select C5apep_W255, C5apep_C2 or C5apep_C3 or C5apep_C4 or C5apep_C5 or C5apep_C7 and chain X and resi 222 and name CA
show spheres, C5apep_W255
set sphere_scale, 0.68, C5apep_W255
color hotpink, C5apep_W255
label C5apep_W255, "W255"
select C5apep_Y258, C5apep_C2 or C5apep_C3 or C5apep_C4 or C5apep_C5 or C5apep_C7 and chain X and resi 225 and name CA
show spheres, C5apep_Y258
set sphere_scale, 0.65, C5apep_Y258
color magenta, C5apep_Y258
label C5apep_Y258, "Y258"
zoom
orient
