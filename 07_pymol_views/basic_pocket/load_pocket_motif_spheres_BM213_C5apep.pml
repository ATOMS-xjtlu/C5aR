reinitialize
cd E:/work/modeling/c5ar/md/article/structure/01_pocket_closeup

# Motif-only sphere view for BM213 rep1 2000 ns vs C5apep rep3 2000 ns.
# Residue mapping used here:
# BM213: M120=87, R134=101, Y222=189, F251=218, W255=222, N296=263, Y300=267
# C5apep: M120=89, R134=103, Y222=191, F251=220, W255=224, N296=265, Y300=269

set_color c5r_cyan, [0.00, 0.72, 0.82]
set_color c5r_orange, [0.92, 0.44, 0.08]

load ../../../gpcr/bm213/bm2000ns.pdb, BM213_rep1
load ../../../gpcr/c5apep/third/2000ns.pdb, C5apep_rep3

# Keep receptor only.
remove BM213_rep1 and not polymer.protein
remove C5apep_rep3 and not polymer.protein

align C5apep_rep3 and polymer.protein, BM213_rep1 and polymer.protein

hide everything

select BM213_motif_spheres, BM213_rep1 and name CA and resi 87+101+189+218+222+263+267
select C5apep_motif_spheres, C5apep_rep3 and name CA and resi 89+103+191+220+224+265+269
select motif_spheres, BM213_motif_spheres or C5apep_motif_spheres

show spheres, motif_spheres
set sphere_scale, 0.75
color c5r_cyan, BM213_motif_spheres
color c5r_orange, C5apep_motif_spheres

set label_size, 14
set label_color, black
set label_outline_color, white

label (BM213_rep1 and name CA and resi 87), "BM213 M120 (zipper)"
label (BM213_rep1 and name CA and resi 101), "BM213 R134 (DRY)"
label (BM213_rep1 and name CA and resi 189), "BM213 Y222 (PIF)"
label (BM213_rep1 and name CA and resi 218), "BM213 F251 (PIF)"
label (BM213_rep1 and name CA and resi 222), "BM213 W255 (CWxP)"
label (BM213_rep1 and name CA and resi 263), "BM213 N296 (NPxxY)"
label (BM213_rep1 and name CA and resi 267), "BM213 Y300 (NPxxY)"

label (C5apep_rep3 and name CA and resi 89), "C5apep M120 (zipper)"
label (C5apep_rep3 and name CA and resi 103), "C5apep R134 (DRY)"
label (C5apep_rep3 and name CA and resi 191), "C5apep Y222 (PIF)"
label (C5apep_rep3 and name CA and resi 220), "C5apep F251 (PIF)"
label (C5apep_rep3 and name CA and resi 224), "C5apep W255 (CWxP)"
label (C5apep_rep3 and name CA and resi 265), "C5apep N296 (NPxxY)"
label (C5apep_rep3 and name CA and resi 269), "C5apep Y300 (NPxxY)"

set bg_rgb, white
set ray_opaque_background, off
orient motif_spheres
zoom motif_spheres, 11
