# C5apep_vs_apo GetContacts node/edge delta
# Open in Windows PyMOL with File > Run Script
reinitialize
load C5apep_md-initial.pdb, receptor
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
set_color delta_positive, [0.10, 0.70, 0.72]
set_color delta_negative, [0.95, 0.55, 0.12]
set_color delta_neutral, [0.95, 0.95, 0.95]
set_color edge_positive, [0.00, 0.55, 0.60]
set_color edge_negative, [0.90, 0.42, 0.00]
set_color key_gray, [0.25, 0.25, 0.25]
set label_size, 16
set label_color, black
set dash_gap, 0.20
set dash_radius, 0.07
set dash_round_ends, on

# Positive color: C5apep higher
# Negative color: apo higher

# Key residue nodes. Cyan means positive state higher. Orange means comparator higher.
select all_key_nodes, none
select node_I91_60, receptor and resi 60
show sticks, node_I91_60
color delta_positive, node_I91_60
show spheres, node_I91_60 and name CA
set sphere_scale, 0.481, node_I91_60 and name CA
label node_I91_60 and name CA, "I91"
select all_key_nodes, all_key_nodes or node_I91_60
select node_W102_71, receptor and resi 71
show sticks, node_W102_71
color delta_negative, node_W102_71
show spheres, node_W102_71 and name CA
set sphere_scale, 0.420, node_W102_71 and name CA
label node_W102_71 and name CA, "W102"
select all_key_nodes, all_key_nodes or node_W102_71
select node_I116_85, receptor and resi 85
show sticks, node_I116_85
color delta_positive, node_I116_85
show spheres, node_I116_85 and name CA
set sphere_scale, 0.461, node_I116_85 and name CA
label node_I116_85 and name CA, "I116"
select all_key_nodes, all_key_nodes or node_I116_85
select node_M120_89, receptor and resi 89
show sticks, node_M120_89
color delta_negative, node_M120_89
show spheres, node_M120_89 and name CA
set sphere_scale, 0.438, node_M120_89 and name CA
label node_M120_89 and name CA, "M120"
select all_key_nodes, all_key_nodes or node_M120_89
select node_R134_103, receptor and resi 103
show sticks, node_R134_103
color delta_negative, node_R134_103
show spheres, node_R134_103 and name CA
set sphere_scale, 0.334, node_R134_103 and name CA
label node_R134_103 and name CA, "R134"
select all_key_nodes, all_key_nodes or node_R134_103
select node_S171_140, receptor and resi 140
show sticks, node_S171_140
color delta_neutral, node_S171_140
show spheres, node_S171_140 and name CA
set sphere_scale, 0.250, node_S171_140 and name CA
label node_S171_140 and name CA, "S171"
select all_key_nodes, all_key_nodes or node_S171_140
select node_Y222_191, receptor and resi 191
show sticks, node_Y222_191
color delta_negative, node_Y222_191
show spheres, node_Y222_191 and name CA
set sphere_scale, 0.457, node_Y222_191 and name CA
label node_Y222_191 and name CA, "Y222"
select all_key_nodes, all_key_nodes or node_Y222_191
select node_F251_220, receptor and resi 220
show sticks, node_F251_220
color delta_negative, node_F251_220
show spheres, node_F251_220 and name CA
set sphere_scale, 0.324, node_F251_220 and name CA
label node_F251_220 and name CA, "F251"
select all_key_nodes, all_key_nodes or node_F251_220
select node_W255_224, receptor and resi 224
show sticks, node_W255_224
color delta_positive, node_W255_224
show spheres, node_W255_224 and name CA
set sphere_scale, 0.485, node_W255_224 and name CA
label node_W255_224 and name CA, "W255"
select all_key_nodes, all_key_nodes or node_W255_224
select node_Y258_227, receptor and resi 227
show sticks, node_Y258_227
color delta_positive, node_Y258_227
show spheres, node_Y258_227 and name CA
set sphere_scale, 0.257, node_Y258_227 and name CA
label node_Y258_227 and name CA, "Y258"
select all_key_nodes, all_key_nodes or node_Y258_227
select node_Y290_259, receptor and resi 259
show sticks, node_Y290_259
color delta_negative, node_Y290_259
show spheres, node_Y290_259 and name CA
set sphere_scale, 0.393, node_Y290_259 and name CA
label node_Y290_259 and name CA, "Y290"
select all_key_nodes, all_key_nodes or node_Y290_259
select node_N292_261, receptor and resi 261
show sticks, node_N292_261
color delta_neutral, node_N292_261
show spheres, node_N292_261 and name CA
set sphere_scale, 0.250, node_N292_261 and name CA
label node_N292_261 and name CA, "N292"
select all_key_nodes, all_key_nodes or node_N292_261
select node_Y300_269, receptor and resi 269
show sticks, node_Y300_269
color delta_positive, node_Y300_269
show spheres, node_Y300_269 and name CA
set sphere_scale, 0.317, node_Y300_269 and name CA
label node_Y300_269 and name CA, "Y300"
select all_key_nodes, all_key_nodes or node_Y300_269

orient all_key_nodes
zoom all_key_nodes, 12

# Top delta edges shown as CA-CA dashed links.
select all_delta_edges, none
distance edge_1_I97_Y222, receptor and resi 99 and name CA, receptor and resi 191 and name CA
hide labels, edge_1_I97_Y222
color edge_negative, edge_1_I97_Y222
set dash_width, 5.66, edge_1_I97_Y222
select all_delta_edges, all_delta_edges or receptor and resi 99+191
distance edge_2_I97_V211, receptor and resi 99 and name CA, receptor and resi 213 and name CA
hide labels, edge_2_I97_V211
color edge_negative, edge_2_I97_V211
set dash_width, 5.60, edge_2_I97_V211
select all_delta_edges, all_delta_edges or receptor and resi 99+213
distance edge_3_R145_Y159, receptor and resi 147 and name CA, receptor and resi 161 and name CA
hide labels, edge_3_R145_Y159
color edge_positive, edge_3_R145_Y159
set dash_width, 4.38, edge_3_R145_Y159
select all_delta_edges, all_delta_edges or receptor and resi 147+161
distance edge_4_R134_V211, receptor and resi 103 and name CA, receptor and resi 213 and name CA
hide labels, edge_4_R134_V211
color edge_negative, edge_4_R134_V211
set dash_width, 4.11, edge_4_R134_V211
select all_delta_edges, all_delta_edges or receptor and resi 103+213
distance edge_5_I91_Q226, receptor and resi 93 and name CA, receptor and resi 228 and name CA
hide labels, edge_5_I91_Q226
color edge_positive, edge_5_I91_Q226
set dash_width, 4.00, edge_5_I91_Q226
select all_delta_edges, all_delta_edges or receptor and resi 93+228
distance edge_6_I91_I116, receptor and resi 60 and name CA, receptor and resi 85 and name CA
hide labels, edge_6_I91_I116
color edge_positive, edge_6_I91_I116
set dash_width, 3.64, edge_6_I91_I116
select all_delta_edges, all_delta_edges or receptor and resi 60+85
distance edge_7_F251_Q226, receptor and resi 220 and name CA, receptor and resi 228 and name CA
hide labels, edge_7_F251_Q226
color edge_negative, edge_7_F251_Q226
set dash_width, 3.21, edge_7_F251_Q226
select all_delta_edges, all_delta_edges or receptor and resi 220+228
distance edge_8_R134_T207, receptor and resi 103 and name CA, receptor and resi 209 and name CA
hide labels, edge_8_R134_T207
color edge_negative, edge_8_R134_T207
set dash_width, 3.19, edge_8_R134_T207
select all_delta_edges, all_delta_edges or receptor and resi 103+209
distance edge_9_I116_Y258, receptor and resi 85 and name CA, receptor and resi 227 and name CA
hide labels, edge_9_I116_Y258
color edge_positive, edge_9_I116_Y258
set dash_width, 3.16, edge_9_I116_Y258
select all_delta_edges, all_delta_edges or receptor and resi 85+227
distance edge_10_W102_L79, receptor and resi 71 and name CA, receptor and resi 81 and name CA
hide labels, edge_10_W102_L79
color edge_negative, edge_10_W102_L79
set dash_width, 3.07, edge_10_W102_L79
select all_delta_edges, all_delta_edges or receptor and resi 71+81
distance edge_11_L59_L79, receptor and resi 61 and name CA, receptor and resi 81 and name CA
hide labels, edge_11_L59_L79
color edge_negative, edge_11_L59_L79
set dash_width, 3.06, edge_11_L59_L79
select all_delta_edges, all_delta_edges or receptor and resi 61+81
distance edge_12_Y88_W255, receptor and resi 90 and name CA, receptor and resi 224 and name CA
hide labels, edge_12_Y88_W255
color edge_positive, edge_12_Y88_W255
set dash_width, 2.97, edge_12_Y88_W255
select all_delta_edges, all_delta_edges or receptor and resi 90+224
distance edge_13_F42_A270, receptor and resi 44 and name CA, receptor and resi 272 and name CA
hide labels, edge_13_F42_A270
color edge_negative, edge_13_F42_A270
set dash_width, 2.90, edge_13_F42_A270
select all_delta_edges, all_delta_edges or receptor and resi 44+272
distance edge_14_H161_E236, receptor and resi 163 and name CA, receptor and resi 238 and name CA
hide labels, edge_14_H161_E236
color edge_positive, edge_14_H161_E236
set dash_width, 2.87, edge_14_H161_E236
select all_delta_edges, all_delta_edges or receptor and resi 163+238
distance edge_15_R145_Y148, receptor and resi 147 and name CA, receptor and resi 150 and name CA
hide labels, edge_15_R145_Y148
color edge_positive, edge_15_R145_Y148
set dash_width, 2.84, edge_15_R145_Y148
select all_delta_edges, all_delta_edges or receptor and resi 147+150


set stick_radius, 0.16, all_key_nodes
show cartoon, receptor
zoom all_key_nodes, 10
# Optional rendering commands:
# ray 2400, 1800
# png C5apep_vs_apo_GetContacts_network.png, dpi=300

