reinitialize
cd E:/work/modeling/c5ar/md/article/structure/01_pocket_closeup

set_color c5r_gray, [0.80, 0.80, 0.80]
set_color c5r_slate, [0.37, 0.45, 0.58]
set_color c5r_blue, [0.10, 0.30, 0.70]
set_color c5r_teal, [0.00, 0.45, 0.50]
set_color c5r_wheat, [0.86, 0.76, 0.55]

load BM213_pocket_initial.pdb, BM213_initial
load BM213_pocket_rep2_2000ns.pdb, BM213_rep2
load BM213_pocket_rep3_2000ns.pdb, BM213_rep3

hide everything
show cartoon, BM213_initial or BM213_rep2 or BM213_rep3
set cartoon_transparency, 0.55, BM213_initial
set cartoon_transparency, 0.30, BM213_rep2
set cartoon_transparency, 0.10, BM213_rep3
color c5r_gray, BM213_initial
color c5r_slate, BM213_rep2
color c5r_blue, BM213_rep3

select BM213_iwi, (BM213_initial or BM213_rep2 or BM213_rep3) and resi 58+69+83
select BM213_switches, (BM213_initial or BM213_rep2 or BM213_rep3) and resi 87+101+138+189+222+225+257+259+263+267
select BM213_ligand, (BM213_initial or BM213_rep2 or BM213_rep3) and resi 282-289
select BM213_closeup, BM213_iwi or BM213_switches or BM213_ligand

show sticks, BM213_closeup
color c5r_wheat, BM213_iwi
color c5r_teal, BM213_switches
color orange, BM213_ligand

set stick_radius, 0.18
set valence, 0
hide everything, hydrogens
set label_size, 16
set label_color, black
set label_outline_color, white

label (BM213_rep3 and name CA and resi 87), "M120^3.36 / resi 87"
label (BM213_rep3 and name CA and resi 101), "R134^3.50 / resi 101"
label (BM213_rep3 and name CA and resi 138), "S171^4.60 / resi 138"
label (BM213_rep3 and name CA and resi 189), "Y222^5.58 / resi 189"
label (BM213_rep3 and name CA and resi 222), "W255^6.48 / resi 222"
label (BM213_rep3 and name CA and resi 225), "Y258^6.51 / resi 225"
label (BM213_rep3 and name CA and resi 257), "Y290^7.43 / resi 257"
label (BM213_rep3 and name CA and resi 259), "N292^7.45 / resi 259"
label (BM213_rep3 and name CA and resi 263), "N296^7.49 / resi 263"
label (BM213_rep3 and name CA and resi 267), "Y300^7.53 / resi 267"

set bg_rgb, white
set ray_opaque_background, off
orient BM213_closeup
zoom BM213_closeup, 7
