# C5apep candidate activation-like communication routes
# Candidate communication routes inferred from AlloViz-derived networks.
# Open in Windows PyMOL with File > Run Script.
reinitialize
load "E:/work/modeling/c5ar/md/article/analysis/alloviz/figures/pymol/C5apep_md-initial.pdb", receptor
remove receptor and not polymer.protein
hide everything
bg_color white
set ray_opaque_background, off
set antialias, 2
set ambient, 0.35
set cartoon_fancy_helices, on
show cartoon, receptor
color gray85, receptor
set cartoon_transparency, 0.20, receptor
set_color route_main, [0.95, 0.55, 0.12]
set_color route_intermediate, [1.00, 0.78, 0.48]
set_color route_supplement, [0.70, 0.35, 0.00]
set_color route_label, [0.05, 0.05, 0.05]
set label_color, route_label
set label_size, 22
set dash_gap, 0.00
set dash_radius, 0.10
set dash_round_ends, on
# Line style controls:
# - Each route segment below uses edge-score-derived dash_width/dash_radius.
# - Larger dash_width/dash_radius means stronger mean edge score.
# - set dash_gap, 0.00 gives a solid-looking line.

# Route residue nodes.
select active_route_nodes, none
select node_F75_44, receptor and resi 44
show spheres, node_F75_44 and name CA
set sphere_scale, 0.48, node_F75_44 and name CA
color route_intermediate, node_F75_44 and name CA
label node_F75_44 and name CA, "F75"
select active_route_nodes, active_route_nodes or node_F75_44
select node_I91_60, receptor and resi 60
show spheres, node_I91_60 and name CA
set sphere_scale, 0.72, node_I91_60 and name CA
color route_supplement, node_I91_60 and name CA
label node_I91_60 and name CA, "I91"
select active_route_nodes, active_route_nodes or node_I91_60
select node_I116_85, receptor and resi 85
show spheres, node_I116_85 and name CA
set sphere_scale, 0.48, node_I116_85 and name CA
color route_supplement, node_I116_85 and name CA
label node_I116_85 and name CA, "I116"
select active_route_nodes, active_route_nodes or node_I116_85
select node_M120_89, receptor and resi 89
show spheres, node_M120_89 and name CA
set sphere_scale, 0.72, node_M120_89 and name CA
color route_supplement, node_M120_89 and name CA
label node_M120_89 and name CA, "M120"
select active_route_nodes, active_route_nodes or node_M120_89
select node_S123_92, receptor and resi 92
show spheres, node_S123_92 and name CA
set sphere_scale, 0.48, node_S123_92 and name CA
color route_intermediate, node_S123_92 and name CA
label node_S123_92 and name CA, "S123"
select active_route_nodes, active_route_nodes or node_S123_92
select node_R134_103, receptor and resi 103
show spheres, node_R134_103 and name CA
set sphere_scale, 0.72, node_R134_103 and name CA
color route_supplement, node_R134_103 and name CA
label node_R134_103 and name CA, "R134"
select active_route_nodes, active_route_nodes or node_R134_103
select node_S171_140, receptor and resi 140
show spheres, node_S171_140 and name CA
set sphere_scale, 0.72, node_S171_140 and name CA
color route_main, node_S171_140 and name CA
label node_S171_140 and name CA, "S171"
select active_route_nodes, active_route_nodes or node_S171_140
select node_R206_175, receptor and resi 175
show spheres, node_R206_175 and name CA
set sphere_scale, 0.72, node_R206_175 and name CA
color route_main, node_R206_175 and name CA
label node_R206_175 and name CA, "R206"
select active_route_nodes, active_route_nodes or node_R206_175
select node_Y222_191, receptor and resi 191
show spheres, node_Y222_191 and name CA
set sphere_scale, 0.72, node_Y222_191 and name CA
color route_main, node_Y222_191 and name CA
label node_Y222_191 and name CA, "Y222"
select active_route_nodes, active_route_nodes or node_Y222_191
select node_T229_198, receptor and resi 198
show spheres, node_T229_198 and name CA
set sphere_scale, 0.48, node_T229_198 and name CA
color route_supplement, node_T229_198 and name CA
label node_T229_198 and name CA, "T229"
select active_route_nodes, active_route_nodes or node_T229_198
select node_F251_220, receptor and resi 220
show spheres, node_F251_220 and name CA
set sphere_scale, 0.72, node_F251_220 and name CA
color route_main, node_F251_220 and name CA
label node_F251_220 and name CA, "F251"
select active_route_nodes, active_route_nodes or node_F251_220
select node_W255_224, receptor and resi 224
show spheres, node_W255_224 and name CA
set sphere_scale, 0.72, node_W255_224 and name CA
color route_main, node_W255_224 and name CA
label node_W255_224 and name CA, "W255"
select active_route_nodes, active_route_nodes or node_W255_224
select node_Y258_227, receptor and resi 227
show spheres, node_Y258_227 and name CA
set sphere_scale, 0.48, node_Y258_227 and name CA
color route_supplement, node_Y258_227 and name CA
label node_Y258_227 and name CA, "Y258"
select active_route_nodes, active_route_nodes or node_Y258_227
select node_N296_265, receptor and resi 265
show spheres, node_N296_265 and name CA
set sphere_scale, 0.72, node_N296_265 and name CA
color route_main, node_N296_265 and name CA
label node_N296_265 and name CA, "N296"
select active_route_nodes, active_route_nodes or node_N296_265
select node_Y300_269, receptor and resi 269
show spheres, node_Y300_269 and name CA
set sphere_scale, 0.72, node_Y300_269 and name CA
color route_main, node_Y300_269 and name CA
label node_Y300_269 and name CA, "Y300"
select active_route_nodes, active_route_nodes or node_Y300_269
select node_R310_279, receptor and resi 279
show spheres, node_R310_279 and name CA
set sphere_scale, 0.48, node_R310_279 and name CA
color route_intermediate, node_R310_279 and name CA
label node_R310_279 and name CA, "R310"
select active_route_nodes, active_route_nodes or node_R310_279

# Route segments as thick solid-looking CA-CA links.
# PIF/CWxP local route: F251 -> S123 -> W255
# segment F251-S123; method=GetContacts; mean_score=0.0441607; support_repeats=2
distance route_1_1_F251_W255, receptor and resi 220 and name CA, receptor and resi 92 and name CA
hide labels, route_1_1_F251_W255
color route_main, route_1_1_F251_W255
set dash_width, 3.57, route_1_1_F251_W255
set dash_radius, 0.07, route_1_1_F251_W255
set dash_gap, 0.00, route_1_1_F251_W255
# segment S123-W255; method=GetContacts; mean_score=0.0337432; support_repeats=3
distance route_1_2_F251_W255, receptor and resi 92 and name CA, receptor and resi 224 and name CA
hide labels, route_1_2_F251_W255
color route_main, route_1_2_F251_W255
set dash_width, 3.37, route_1_2_F251_W255
set dash_radius, 0.07, route_1_2_F251_W255
set dash_gap, 0.00, route_1_2_F251_W255
# TM7/NPxxY-adjacent contact route: N296 -> F75 -> R310 -> Y300
# segment N296-F75; method=GetContacts; mean_score=0.162148; support_repeats=3
distance route_2_1_N296_Y300, receptor and resi 265 and name CA, receptor and resi 44 and name CA
hide labels, route_2_1_N296_Y300
color route_main, route_2_1_N296_Y300
set dash_width, 5.79, route_2_1_N296_Y300
set dash_radius, 0.12, route_2_1_N296_Y300
set dash_gap, 0.00, route_2_1_N296_Y300
# segment F75-R310; method=GetContacts; mean_score=0.173238; support_repeats=3
distance route_2_2_N296_Y300, receptor and resi 44 and name CA, receptor and resi 279 and name CA
hide labels, route_2_2_N296_Y300
color route_main, route_2_2_N296_Y300
set dash_width, 6.00, route_2_2_N296_Y300
set dash_radius, 0.12, route_2_2_N296_Y300
set dash_gap, 0.00, route_2_2_N296_Y300
# segment R310-Y300; method=GetContacts; mean_score=0.102024; support_repeats=3
distance route_2_3_N296_Y300, receptor and resi 279 and name CA, receptor and resi 269 and name CA
hide labels, route_2_3_N296_Y300
color route_main, route_2_3_N296_Y300
set dash_width, 4.66, route_2_3_N296_Y300
set dash_radius, 0.09, route_2_3_N296_Y300
set dash_gap, 0.00, route_2_3_N296_Y300
# IWI-to-TM3/CWxP contact route: I91 -> I116 -> Y258 -> M120
# segment I91-I116; method=GetContacts; mean_score=0.205224; support_repeats=3
distance route_3_1_I91_M120, receptor and resi 60 and name CA, receptor and resi 85 and name CA
hide labels, route_3_1_I91_M120
color route_supplement, route_3_1_I91_M120
set dash_width, 6.60, route_3_1_I91_M120
set dash_radius, 0.13, route_3_1_I91_M120
set dash_gap, 0.00, route_3_1_I91_M120
# segment I116-Y258; method=GetContacts; mean_score=0.333005; support_repeats=2
distance route_3_2_I91_M120, receptor and resi 85 and name CA, receptor and resi 227 and name CA
hide labels, route_3_2_I91_M120
color route_supplement, route_3_2_I91_M120
set dash_width, 9.00, route_3_2_I91_M120
set dash_radius, 0.18, route_3_2_I91_M120
set dash_gap, 0.00, route_3_2_I91_M120
# segment Y258-M120; method=GetContacts; mean_score=0.177545; support_repeats=3
distance route_3_3_I91_M120, receptor and resi 227 and name CA, receptor and resi 89 and name CA
hide labels, route_3_3_I91_M120
color route_supplement, route_3_3_I91_M120
set dash_width, 6.08, route_3_3_I91_M120
set dash_radius, 0.12, route_3_3_I91_M120
set dash_gap, 0.00, route_3_3_I91_M120
# TM3-to-NPxxY correlation route: M120 -> Y300
# segment M120-Y300; method=correlationplus_CA_Pear; mean_score=8.35359e-05; support_repeats=3
distance route_4_1_M120_Y300, receptor and resi 89 and name CA, receptor and resi 269 and name CA
hide labels, route_4_1_M120_Y300
color route_supplement, route_4_1_M120_Y300
set dash_width, 3.09, route_4_1_M120_Y300
set dash_radius, 0.06, route_4_1_M120_Y300
set dash_gap, 0.00, route_4_1_M120_Y300
# DRY-like-to-TM5 correlation route: R134 -> T229 -> Y222
# segment R134-T229; method=correlationplus_CA_Pear; mean_score=5.84751e-05; support_repeats=3
distance route_5_1_R134_Y222, receptor and resi 103 and name CA, receptor and resi 198 and name CA
hide labels, route_5_1_R134_Y222
color route_supplement, route_5_1_R134_Y222
set dash_width, 3.00, route_5_1_R134_Y222
set dash_radius, 0.06, route_5_1_R134_Y222
set dash_gap, 0.00, route_5_1_R134_Y222
# segment T229-Y222; method=correlationplus_CA_Pear; mean_score=0.00167907; support_repeats=3
distance route_5_2_R134_Y222, receptor and resi 198 and name CA, receptor and resi 191 and name CA
hide labels, route_5_2_R134_Y222
color route_supplement, route_5_2_R134_Y222
set dash_width, 9.00, route_5_2_R134_Y222
set dash_radius, 0.18, route_5_2_R134_Y222
set dash_gap, 0.00, route_5_2_R134_Y222
# peptide-contact route: S171 -> R206
# segment S171-R206; method=GetContacts; mean_score=0.0138473; support_repeats=3
distance route_6_1_S171_R206, receptor and resi 140 and name CA, receptor and resi 175 and name CA
hide labels, route_6_1_S171_R206
color route_main, route_6_1_S171_R206
set dash_width, 3.00, route_6_1_S171_R206
set dash_radius, 0.06, route_6_1_S171_R206
set dash_gap, 0.00, route_6_1_S171_R206
# Y222-F251 link: Y222 -> F251
# segment Y222-F251; method=GetContacts; mean_score=0.323947; support_repeats=1
distance route_7_1_Y222_F251, receptor and resi 191 and name CA, receptor and resi 220 and name CA
hide labels, route_7_1_Y222_F251
color route_main, route_7_1_Y222_F251
set dash_width, 8.83, route_7_1_Y222_F251
set dash_radius, 0.18, route_7_1_Y222_F251
set dash_gap, 0.00, route_7_1_Y222_F251

orient active_route_nodes
zoom active_route_nodes, 12
# Optional rendering commands:
# ray 2400, 1800
# png C5apep_activation_route.png, dpi=300

