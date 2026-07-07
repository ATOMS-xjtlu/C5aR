# BM213_vs_apo GetContacts node/edge delta
# Open in Windows PyMOL with File > Run Script
reinitialize
load BM213_md-initial.pdb, receptor
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

# Positive color: BM213 higher
# Negative color: apo higher

# Key residue nodes. Cyan means positive state higher. Orange means comparator higher.
select all_key_nodes, none
select node_I91_58, receptor and resi 58
show sticks, node_I91_58
color delta_positive, node_I91_58
show spheres, node_I91_58 and name CA
set sphere_scale, 0.438, node_I91_58 and name CA
label node_I91_58 and name CA, "I91"
select all_key_nodes, all_key_nodes or node_I91_58
select node_W102_69, receptor and resi 69
show sticks, node_W102_69
color delta_positive, node_W102_69
show spheres, node_W102_69 and name CA
set sphere_scale, 0.394, node_W102_69 and name CA
label node_W102_69 and name CA, "W102"
select all_key_nodes, all_key_nodes or node_W102_69
select node_I116_83, receptor and resi 83
show sticks, node_I116_83
color delta_positive, node_I116_83
show spheres, node_I116_83 and name CA
set sphere_scale, 0.289, node_I116_83 and name CA
label node_I116_83 and name CA, "I116"
select all_key_nodes, all_key_nodes or node_I116_83
select node_M120_87, receptor and resi 87
show sticks, node_M120_87
color delta_positive, node_M120_87
show spheres, node_M120_87 and name CA
set sphere_scale, 0.310, node_M120_87 and name CA
label node_M120_87 and name CA, "M120"
select all_key_nodes, all_key_nodes or node_M120_87
select node_R134_101, receptor and resi 101
show sticks, node_R134_101
color delta_negative, node_R134_101
show spheres, node_R134_101 and name CA
set sphere_scale, 0.452, node_R134_101 and name CA
label node_R134_101 and name CA, "R134"
select all_key_nodes, all_key_nodes or node_R134_101
select node_S171_138, receptor and resi 138
show sticks, node_S171_138
color delta_positive, node_S171_138
show spheres, node_S171_138 and name CA
set sphere_scale, 0.320, node_S171_138 and name CA
label node_S171_138 and name CA, "S171"
select all_key_nodes, all_key_nodes or node_S171_138
select node_Y222_189, receptor and resi 189
show sticks, node_Y222_189
color delta_negative, node_Y222_189
show spheres, node_Y222_189 and name CA
set sphere_scale, 0.427, node_Y222_189 and name CA
label node_Y222_189 and name CA, "Y222"
select all_key_nodes, all_key_nodes or node_Y222_189
select node_F251_218, receptor and resi 218
show sticks, node_F251_218
color delta_negative, node_F251_218
show spheres, node_F251_218 and name CA
set sphere_scale, 0.339, node_F251_218 and name CA
label node_F251_218 and name CA, "F251"
select all_key_nodes, all_key_nodes or node_F251_218
select node_W255_222, receptor and resi 222
show sticks, node_W255_222
color delta_positive, node_W255_222
show spheres, node_W255_222 and name CA
set sphere_scale, 0.614, node_W255_222 and name CA
label node_W255_222 and name CA, "W255"
select all_key_nodes, all_key_nodes or node_W255_222
select node_Y258_225, receptor and resi 225
show sticks, node_Y258_225
color delta_positive, node_Y258_225
show spheres, node_Y258_225 and name CA
set sphere_scale, 0.352, node_Y258_225 and name CA
label node_Y258_225 and name CA, "Y258"
select all_key_nodes, all_key_nodes or node_Y258_225
select node_Y290_257, receptor and resi 257
show sticks, node_Y290_257
color delta_positive, node_Y290_257
show spheres, node_Y290_257 and name CA
set sphere_scale, 0.423, node_Y290_257 and name CA
label node_Y290_257 and name CA, "Y290"
select all_key_nodes, all_key_nodes or node_Y290_257
select node_N292_259, receptor and resi 259
show sticks, node_N292_259
color delta_neutral, node_N292_259
show spheres, node_N292_259 and name CA
set sphere_scale, 0.250, node_N292_259 and name CA
label node_N292_259 and name CA, "N292"
select all_key_nodes, all_key_nodes or node_N292_259
select node_Y300_267, receptor and resi 267
show sticks, node_Y300_267
color delta_negative, node_Y300_267
show spheres, node_Y300_267 and name CA
set sphere_scale, 0.262, node_Y300_267 and name CA
label node_Y300_267 and name CA, "Y300"
select all_key_nodes, all_key_nodes or node_Y300_267

orient all_key_nodes
zoom all_key_nodes, 12

# Top delta edges shown as CA-CA dashed links.
select all_delta_edges, none
distance edge_1_I97_Y222, receptor and resi 97 and name CA, receptor and resi 189 and name CA
hide labels, edge_1_I97_Y222
color edge_negative, edge_1_I97_Y222
set dash_width, 5.66, edge_1_I97_Y222
select all_delta_edges, all_delta_edges or receptor and resi 97+189
distance edge_2_F251_Q226, receptor and resi 218 and name CA, receptor and resi 226 and name CA
hide labels, edge_2_F251_Q226
color edge_negative, edge_2_F251_Q226
set dash_width, 5.09, edge_2_F251_Q226
select all_delta_edges, all_delta_edges or receptor and resi 218+226
distance edge_3_R134_V211, receptor and resi 101 and name CA, receptor and resi 211 and name CA
hide labels, edge_3_R134_V211
color edge_negative, edge_3_R134_V211
set dash_width, 4.64, edge_3_R134_V211
select all_delta_edges, all_delta_edges or receptor and resi 101+211
distance edge_4_H68_K246, receptor and resi 68 and name CA, receptor and resi 246 and name CA
hide labels, edge_4_H68_K246
color edge_positive, edge_4_H68_K246
set dash_width, 3.88, edge_4_H68_K246
select all_delta_edges, all_delta_edges or receptor and resi 68+246
distance edge_5_M120_W255, receptor and resi 87 and name CA, receptor and resi 222 and name CA
hide labels, edge_5_M120_W255
color edge_positive, edge_5_M120_W255
set dash_width, 3.81, edge_5_M120_W255
select all_delta_edges, all_delta_edges or receptor and resi 87+222
distance edge_6_Y222_F251, receptor and resi 189 and name CA, receptor and resi 218 and name CA
hide labels, edge_6_Y222_F251
color edge_negative, edge_6_Y222_F251
set dash_width, 3.76, edge_6_Y222_F251
select all_delta_edges, all_delta_edges or receptor and resi 189+218
distance edge_7_L59_W102, receptor and resi 59 and name CA, receptor and resi 69 and name CA
hide labels, edge_7_L59_W102
color edge_positive, edge_7_L59_W102
set dash_width, 3.70, edge_7_L59_W102
select all_delta_edges, all_delta_edges or receptor and resi 59+69
distance edge_8_M120_R173, receptor and resi 87 and name CA, receptor and resi 173 and name CA
hide labels, edge_8_M120_R173
color edge_negative, edge_8_M120_R173
set dash_width, 3.52, edge_8_M120_R173
select all_delta_edges, all_delta_edges or receptor and resi 87+173
distance edge_9_M120_Q226, receptor and resi 87 and name CA, receptor and resi 226 and name CA
hide labels, edge_9_M120_Q226
color edge_negative, edge_9_M120_Q226
set dash_width, 3.51, edge_9_M120_Q226
select all_delta_edges, all_delta_edges or receptor and resi 87+226
distance edge_10_F42_R134, receptor and resi 42 and name CA, receptor and resi 101 and name CA
hide labels, edge_10_F42_R134
color edge_negative, edge_10_F42_R134
set dash_width, 3.42, edge_10_F42_R134
select all_delta_edges, all_delta_edges or receptor and resi 42+101
distance edge_11_R134_T207, receptor and resi 101 and name CA, receptor and resi 207 and name CA
hide labels, edge_11_R134_T207
color edge_negative, edge_11_R134_T207
set dash_width, 3.25, edge_11_R134_T207
select all_delta_edges, all_delta_edges or receptor and resi 101+207
distance edge_12_L84_Y290, receptor and resi 84 and name CA, receptor and resi 257 and name CA
hide labels, edge_12_L84_Y290
color edge_positive, edge_12_L84_Y290
set dash_width, 3.21, edge_12_L84_Y290
select all_delta_edges, all_delta_edges or receptor and resi 84+257
distance edge_13_W102_L79, receptor and resi 69 and name CA, receptor and resi 79 and name CA
hide labels, edge_13_W102_L79
color edge_negative, edge_13_W102_L79
set dash_width, 3.17, edge_13_W102_L79
select all_delta_edges, all_delta_edges or receptor and resi 69+79
distance edge_14_I116_Y258, receptor and resi 83 and name CA, receptor and resi 225 and name CA
hide labels, edge_14_I116_Y258
color edge_positive, edge_14_I116_Y258
set dash_width, 3.14, edge_14_I116_Y258
select all_delta_edges, all_delta_edges or receptor and resi 83+225
distance edge_15_H68_Y148, receptor and resi 68 and name CA, receptor and resi 148 and name CA
hide labels, edge_15_H68_Y148
color edge_positive, edge_15_H68_Y148
set dash_width, 3.13, edge_15_H68_Y148
select all_delta_edges, all_delta_edges or receptor and resi 68+148


set stick_radius, 0.16, all_key_nodes
show cartoon, receptor
zoom all_key_nodes, 10
# Optional rendering commands:
# ray 2400, 1800
# png BM213_vs_apo_GetContacts_network.png, dpi=300

