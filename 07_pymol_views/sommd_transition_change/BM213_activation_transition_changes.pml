# BM213 major activation-core transition changes
# Spheres highlight residues with largest transition-weighted changes.
reinitialize
bg_color white
set cartoon_transparency, 0.25
set sphere_quality, 2
set label_size, 18
load /mnt/e/work/modeling/c5ar/md/article/analysis/sommd/pymol/formal/cluster_01_BM213_rep3.pdb, BM213_C1
hide everything, BM213_C1
show cartoon, BM213_C1
color tv_blue, BM213_C1
load /mnt/e/work/modeling/c5ar/md/article/analysis/sommd/pymol/formal/cluster_04_BM213_rep3.pdb, BM213_C4
hide everything, BM213_C4
show cartoon, BM213_C4
color tv_blue, BM213_C4
load /mnt/e/work/modeling/c5ar/md/article/analysis/sommd/pymol/formal/cluster_05_BM213_rep3.pdb, BM213_C5
hide everything, BM213_C5
show cartoon, BM213_C5
color tv_blue, BM213_C5
load /mnt/e/work/modeling/c5ar/md/article/analysis/sommd/pymol/formal/archive_non_dominant_cluster_pdbs/cluster_07_BM213_rep2.pdb, BM213_C7
hide everything, BM213_C7
show cartoon, BM213_C7
color tv_blue, BM213_C7
load /mnt/e/work/modeling/c5ar/md/article/analysis/sommd/pymol/formal/cluster_08_BM213_rep2.pdb, BM213_C8
hide everything, BM213_C8
show cartoon, BM213_C8
color tv_blue, BM213_C8
align BM213_C4 and chain X and name CA, BM213_C1 and chain X and name CA
align BM213_C5 and chain X and name CA, BM213_C1 and chain X and name CA
align BM213_C7 and chain X and name CA, BM213_C1 and chain X and name CA
align BM213_C8 and chain X and name CA, BM213_C1 and chain X and name CA
select BM213_W102, BM213_C1 or BM213_C4 or BM213_C5 or BM213_C7 or BM213_C8 and chain X and resi 69 and name CA
show spheres, BM213_W102
set sphere_scale, 0.71, BM213_W102
color red, BM213_W102
label BM213_W102, "W102"
select BM213_R134, BM213_C1 or BM213_C4 or BM213_C5 or BM213_C7 or BM213_C8 and chain X and resi 101 and name CA
show spheres, BM213_R134
set sphere_scale, 0.64, BM213_R134
color orange, BM213_R134
label BM213_R134, "R134"
select BM213_Y300, BM213_C1 or BM213_C4 or BM213_C5 or BM213_C7 or BM213_C8 and chain X and resi 267 and name CA
show spheres, BM213_Y300
set sphere_scale, 0.63, BM213_Y300
color yelloworange, BM213_Y300
label BM213_Y300, "Y300"
select BM213_F251, BM213_C1 or BM213_C4 or BM213_C5 or BM213_C7 or BM213_C8 and chain X and resi 218 and name CA
show spheres, BM213_F251
set sphere_scale, 0.62, BM213_F251
color tv_yellow, BM213_F251
label BM213_F251, "F251"
select BM213_Y222, BM213_C1 or BM213_C4 or BM213_C5 or BM213_C7 or BM213_C8 and chain X and resi 189 and name CA
show spheres, BM213_Y222
set sphere_scale, 0.61, BM213_Y222
color salmon, BM213_Y222
label BM213_Y222, "Y222"
select BM213_N296, BM213_C1 or BM213_C4 or BM213_C5 or BM213_C7 or BM213_C8 and chain X and resi 263 and name CA
show spheres, BM213_N296
set sphere_scale, 0.59, BM213_N296
color raspberry, BM213_N296
label BM213_N296, "N296"
select BM213_Y258, BM213_C1 or BM213_C4 or BM213_C5 or BM213_C7 or BM213_C8 and chain X and resi 225 and name CA
show spheres, BM213_Y258
set sphere_scale, 0.59, BM213_Y258
color hotpink, BM213_Y258
label BM213_Y258, "Y258"
select BM213_W255, BM213_C1 or BM213_C4 or BM213_C5 or BM213_C7 or BM213_C8 and chain X and resi 222 and name CA
show spheres, BM213_W255
set sphere_scale, 0.58, BM213_W255
color magenta, BM213_W255
label BM213_W255, "W255"
zoom
orient
