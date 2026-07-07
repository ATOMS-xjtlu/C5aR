# BM213 candidate activation-like communication routes
# Candidate communication routes inferred from AlloViz-derived networks.
# Open in Windows PyMOL with File > Run Script.
reinitialize
load "E:/work/modeling/c5ar/md/article/analysis/alloviz/figures/pymol/BM213_md-initial.pdb", receptor
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
set_color route_main, [0.05, 0.65, 0.68]
set_color route_intermediate, [0.55, 0.86, 0.88]
set_color route_supplement, [0.28, 0.40, 0.78]
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
select node_I116_83, receptor and resi 83
show spheres, node_I116_83 and name CA
set sphere_scale, 0.72, node_I116_83 and name CA
color route_main, node_I116_83 and name CA
label node_I116_83 and name CA, "I116"
select active_route_nodes, active_route_nodes or node_I116_83
select node_M120_87, receptor and resi 87
show spheres, node_M120_87 and name CA
set sphere_scale, 0.72, node_M120_87 and name CA
color route_main, node_M120_87 and name CA
label node_M120_87 and name CA, "M120"
select active_route_nodes, active_route_nodes or node_M120_87
select node_L127_94, receptor and resi 94
show spheres, node_L127_94 and name CA
set sphere_scale, 0.48, node_L127_94 and name CA
color route_supplement, node_L127_94 and name CA
label node_L127_94 and name CA, "L127"
select active_route_nodes, active_route_nodes or node_L127_94
select node_R134_101, receptor and resi 101
show spheres, node_R134_101 and name CA
set sphere_scale, 0.72, node_R134_101 and name CA
color route_main, node_R134_101 and name CA
label node_R134_101 and name CA, "R134"
select active_route_nodes, active_route_nodes or node_R134_101
select node_S171_138, receptor and resi 138
show spheres, node_S171_138 and name CA
set sphere_scale, 0.72, node_S171_138 and name CA
color route_main, node_S171_138 and name CA
label node_S171_138 and name CA, "S171"
select active_route_nodes, active_route_nodes or node_S171_138
select node_R206_173, receptor and resi 173
show spheres, node_R206_173 and name CA
set sphere_scale, 0.48, node_R206_173 and name CA
color route_intermediate, node_R206_173 and name CA
label node_R206_173 and name CA, "R206"
select active_route_nodes, active_route_nodes or node_R206_173
select node_Y222_189, receptor and resi 189
show spheres, node_Y222_189 and name CA
set sphere_scale, 0.72, node_Y222_189 and name CA
color route_main, node_Y222_189 and name CA
label node_Y222_189 and name CA, "Y222"
select active_route_nodes, active_route_nodes or node_Y222_189
select node_V243_210, receptor and resi 210
show spheres, node_V243_210 and name CA
set sphere_scale, 0.48, node_V243_210 and name CA
color route_intermediate, node_V243_210 and name CA
label node_V243_210 and name CA, "V243"
select active_route_nodes, active_route_nodes or node_V243_210
select node_F251_218, receptor and resi 218
show spheres, node_F251_218 and name CA
set sphere_scale, 0.72, node_F251_218 and name CA
color route_main, node_F251_218 and name CA
label node_F251_218 and name CA, "F251"
select active_route_nodes, active_route_nodes or node_F251_218
select node_W255_222, receptor and resi 222
show spheres, node_W255_222 and name CA
set sphere_scale, 0.72, node_W255_222 and name CA
color route_main, node_W255_222 and name CA
label node_W255_222 and name CA, "W255"
select active_route_nodes, active_route_nodes or node_W255_222
select node_Y258_225, receptor and resi 225
show spheres, node_Y258_225 and name CA
set sphere_scale, 0.72, node_Y258_225 and name CA
color route_main, node_Y258_225 and name CA
label node_Y258_225 and name CA, "Y258"
select active_route_nodes, active_route_nodes or node_Y258_225
select node_Y290_257, receptor and resi 257
show spheres, node_Y290_257 and name CA
set sphere_scale, 0.72, node_Y290_257 and name CA
color route_main, node_Y290_257 and name CA
label node_Y290_257 and name CA, "Y290"
select active_route_nodes, active_route_nodes or node_Y290_257

# Route segments as thick solid-looking CA-CA links.
# TM3-CWxP direct bridge: M120 -> W255
# segment M120-W255; method=GetContacts; mean_score=0.312767; support_repeats=3
distance route_1_1_M120_W255, receptor and resi 87 and name CA, receptor and resi 222 and name CA
hide labels, route_1_1_M120_W255
color route_main, route_1_1_M120_W255
set dash_width, 7.22, route_1_1_M120_W255
set dash_radius, 0.14, route_1_1_M120_W255
set dash_gap, 0.00, route_1_1_M120_W255
# CWxP microswitch core: W255 -> M120 -> Y258
# segment W255-M120; method=GetContacts; mean_score=0.312767; support_repeats=3
distance route_2_1_W255_Y258, receptor and resi 222 and name CA, receptor and resi 87 and name CA
hide labels, route_2_1_W255_Y258
color route_main, route_2_1_W255_Y258
set dash_width, 7.22, route_2_1_W255_Y258
set dash_radius, 0.14, route_2_1_W255_Y258
set dash_gap, 0.00, route_2_1_W255_Y258
# segment M120-Y258; method=GetContacts; mean_score=0.313167; support_repeats=3
distance route_2_2_W255_Y258, receptor and resi 87 and name CA, receptor and resi 225 and name CA
hide labels, route_2_2_W255_Y258
color route_main, route_2_2_W255_Y258
set dash_width, 7.22, route_2_2_W255_Y258
set dash_radius, 0.14, route_2_2_W255_Y258
set dash_gap, 0.00, route_2_2_W255_Y258
# TM6-to-TM7/IWI route: W255 -> M120 -> Y258 -> I116 -> Y290
# segment W255-M120; method=GetContacts; mean_score=0.312767; support_repeats=3
distance route_3_1_W255_Y290, receptor and resi 222 and name CA, receptor and resi 87 and name CA
hide labels, route_3_1_W255_Y290
color route_main, route_3_1_W255_Y290
set dash_width, 7.22, route_3_1_W255_Y290
set dash_radius, 0.14, route_3_1_W255_Y290
set dash_gap, 0.00, route_3_1_W255_Y290
# segment M120-Y258; method=GetContacts; mean_score=0.313167; support_repeats=3
distance route_3_2_W255_Y290, receptor and resi 87 and name CA, receptor and resi 225 and name CA
hide labels, route_3_2_W255_Y290
color route_main, route_3_2_W255_Y290
set dash_width, 7.22, route_3_2_W255_Y290
set dash_radius, 0.14, route_3_2_W255_Y290
set dash_gap, 0.00, route_3_2_W255_Y290
# segment Y258-I116; method=GetContacts; mean_score=0.330727; support_repeats=2
distance route_3_3_W255_Y290, receptor and resi 225 and name CA, receptor and resi 83 and name CA
hide labels, route_3_3_W255_Y290
color route_main, route_3_3_W255_Y290
set dash_width, 7.47, route_3_3_W255_Y290
set dash_radius, 0.15, route_3_3_W255_Y290
set dash_gap, 0.00, route_3_3_W255_Y290
# segment I116-Y290; method=GetContacts; mean_score=0.22814; support_repeats=3
distance route_3_4_W255_Y290, receptor and resi 83 and name CA, receptor and resi 257 and name CA
hide labels, route_3_4_W255_Y290
color route_main, route_3_4_W255_Y290
set dash_width, 6.01, route_3_4_W255_Y290
set dash_radius, 0.12, route_3_4_W255_Y290
set dash_gap, 0.00, route_3_4_W255_Y290
# DRY-like-to-TM5 route: R134 -> V243 -> Y222
# segment R134-V243; method=GetContacts; mean_score=0.0175911; support_repeats=2
distance route_4_1_R134_Y222, receptor and resi 101 and name CA, receptor and resi 210 and name CA
hide labels, route_4_1_R134_Y222
color route_main, route_4_1_R134_Y222
set dash_width, 3.00, route_4_1_R134_Y222
set dash_radius, 0.06, route_4_1_R134_Y222
set dash_gap, 0.00, route_4_1_R134_Y222
# segment V243-Y222; method=GetContacts; mean_score=0.0543627; support_repeats=2
distance route_4_2_R134_Y222, receptor and resi 210 and name CA, receptor and resi 189 and name CA
hide labels, route_4_2_R134_Y222
color route_main, route_4_2_R134_Y222
set dash_width, 3.53, route_4_2_R134_Y222
set dash_radius, 0.07, route_4_2_R134_Y222
set dash_gap, 0.00, route_4_2_R134_Y222
# pocket-to-CWxP route: S171 -> R206 -> Y258
# segment S171-R206; method=GetContacts; mean_score=0.0543087; support_repeats=3
distance route_5_1_S171_Y258, receptor and resi 138 and name CA, receptor and resi 173 and name CA
hide labels, route_5_1_S171_Y258
color route_main, route_5_1_S171_Y258
set dash_width, 3.52, route_5_1_S171_Y258
set dash_radius, 0.07, route_5_1_S171_Y258
set dash_gap, 0.00, route_5_1_S171_Y258
# segment R206-Y258; method=GetContacts; mean_score=0.0958971; support_repeats=2
distance route_5_2_S171_Y258, receptor and resi 173 and name CA, receptor and resi 225 and name CA
hide labels, route_5_2_S171_Y258
color route_main, route_5_2_S171_Y258
set dash_width, 4.12, route_5_2_S171_Y258
set dash_radius, 0.08, route_5_2_S171_Y258
set dash_gap, 0.00, route_5_2_S171_Y258
# IWI-to-CWxP correlation route: I116 -> Y258
# segment I116-Y258; method=correlationplus_CA_Pear; mean_score=2.54194e-05; support_repeats=3
distance route_6_1_I116_Y258, receptor and resi 83 and name CA, receptor and resi 225 and name CA
hide labels, route_6_1_I116_Y258
color route_supplement, route_6_1_I116_Y258
set dash_width, 3.00, route_6_1_I116_Y258
set dash_radius, 0.06, route_6_1_I116_Y258
set dash_gap, 0.00, route_6_1_I116_Y258
# TM3-to-CWxP correlation route: M120 -> Y258
# segment M120-Y258; method=correlationplus_CA_Pear; mean_score=2.54194e-05; support_repeats=2
distance route_7_1_M120_Y258, receptor and resi 87 and name CA, receptor and resi 225 and name CA
hide labels, route_7_1_M120_Y258
color route_supplement, route_7_1_M120_Y258
set dash_width, 3.00, route_7_1_M120_Y258
set dash_radius, 0.06, route_7_1_M120_Y258
set dash_gap, 0.00, route_7_1_M120_Y258
# TM3-CWxP correlation route: M120 -> W255
# segment M120-W255; method=correlationplus_CA_Pear; mean_score=0.000127097; support_repeats=2
distance route_8_1_M120_W255, receptor and resi 87 and name CA, receptor and resi 222 and name CA
hide labels, route_8_1_M120_W255
color route_supplement, route_8_1_M120_W255
set dash_width, 9.00, route_8_1_M120_W255
set dash_radius, 0.18, route_8_1_M120_W255
set dash_gap, 0.00, route_8_1_M120_W255
# W255-L127-Y222 bridge: W255 -> L127 -> Y222
# segment W255-L127; method=GetContacts; mean_score=0.437672; support_repeats=2
distance route_9_1_W255_L127_Y222, receptor and resi 222 and name CA, receptor and resi 94 and name CA
hide labels, route_9_1_W255_L127_Y222
color route_supplement, route_9_1_W255_L127_Y222
set dash_width, 9.00, route_9_1_W255_L127_Y222
set dash_radius, 0.18, route_9_1_W255_L127_Y222
set dash_gap, 0.00, route_9_1_W255_L127_Y222
# segment L127-Y222; method=GetContacts; mean_score=0.113684; support_repeats=2
distance route_9_2_W255_L127_Y222, receptor and resi 94 and name CA, receptor and resi 189 and name CA
hide labels, route_9_2_W255_L127_Y222
color route_supplement, route_9_2_W255_L127_Y222
set dash_width, 4.37, route_9_2_W255_L127_Y222
set dash_radius, 0.09, route_9_2_W255_L127_Y222
set dash_gap, 0.00, route_9_2_W255_L127_Y222
# W255-L127-R134 bridge: W255 -> L127 -> R134
# segment W255-L127; method=GetContacts; mean_score=0.437672; support_repeats=2
distance route_10_1_W255_L127_R134, receptor and resi 222 and name CA, receptor and resi 94 and name CA
hide labels, route_10_1_W255_L127_R134
color route_supplement, route_10_1_W255_L127_R134
set dash_width, 9.00, route_10_1_W255_L127_R134
set dash_radius, 0.18, route_10_1_W255_L127_R134
set dash_gap, 0.00, route_10_1_W255_L127_R134
# segment L127-R134; method=GetContacts; mean_score=0.308382; support_repeats=1
distance route_10_2_W255_L127_R134, receptor and resi 94 and name CA, receptor and resi 101 and name CA
hide labels, route_10_2_W255_L127_R134
color route_supplement, route_10_2_W255_L127_R134
set dash_width, 7.15, route_10_2_W255_L127_R134
set dash_radius, 0.14, route_10_2_W255_L127_R134
set dash_gap, 0.00, route_10_2_W255_L127_R134
# Y222-F251 link: Y222 -> F251
# segment Y222-F251; method=GetContacts; mean_score=0.0392857; support_repeats=1
distance route_11_1_Y222_F251, receptor and resi 189 and name CA, receptor and resi 218 and name CA
hide labels, route_11_1_Y222_F251
color route_main, route_11_1_Y222_F251
set dash_width, 3.31, route_11_1_Y222_F251
set dash_radius, 0.07, route_11_1_Y222_F251
set dash_gap, 0.00, route_11_1_Y222_F251

orient active_route_nodes
zoom active_route_nodes, 12
# Optional rendering commands:
# ray 2400, 1800
# png BM213_activation_route.png, dpi=300

