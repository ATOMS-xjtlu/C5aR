reinitialize
cd E:/work/modeling/c5ar/md/article/structure/01_pocket_closeup

set_color c5r_gray, [0.80, 0.80, 0.80]
set_color c5r_orange, [0.90, 0.45, 0.10]
set_color c5r_red, [0.65, 0.16, 0.16]
set_color c5r_wheat, [0.86, 0.76, 0.55]

load C5apep_pocket_initial.pdb, C5apep_initial
load C5apep_pocket_rep3_2000ns.pdb, C5apep_rep3

hide everything
show cartoon, C5apep_initial or C5apep_rep3
set cartoon_transparency, 0.55, C5apep_initial
set cartoon_transparency, 0.10, C5apep_rep3
color c5r_gray, C5apep_initial
color c5r_orange, C5apep_rep3

select C5apep_iwi, (C5apep_initial or C5apep_rep3) and resi 60+71+85
select C5apep_switches, (C5apep_initial or C5apep_rep3) and resi 89+103+140+191+224+227+259+261+265+269
select C5apep_ligand, (C5apep_initial or C5apep_rep3) and resi 284
select C5apep_closeup, C5apep_iwi or C5apep_switches or C5apep_ligand

show sticks, C5apep_closeup
color c5r_wheat, C5apep_iwi
color c5r_red, C5apep_switches
color orange, C5apep_ligand

set stick_radius, 0.18
set valence, 0
hide everything, hydrogens
set label_size, 16
set label_color, black
set label_outline_color, white

label (C5apep_rep3 and name CA and resi 89), "M120^3.36 / resi 89"
label (C5apep_rep3 and name CA and resi 103), "R134^3.50 / resi 103"
label (C5apep_rep3 and name CA and resi 140), "S171^4.60 / resi 140"
label (C5apep_rep3 and name CA and resi 191), "Y222^5.58 / resi 191"
label (C5apep_rep3 and name CA and resi 224), "W255^6.48 / resi 224"
label (C5apep_rep3 and name CA and resi 227), "Y258^6.51 / resi 227"
label (C5apep_rep3 and name CA and resi 259), "Y290^7.43 / resi 259"
label (C5apep_rep3 and name CA and resi 261), "N292^7.45 / resi 261"
label (C5apep_rep3 and name CA and resi 265), "N296^7.49 / resi 265"
label (C5apep_rep3 and name CA and resi 269), "Y300^7.53 / resi 269"

set bg_rgb, white
set ray_opaque_background, off
orient C5apep_closeup
zoom C5apep_closeup, 7
