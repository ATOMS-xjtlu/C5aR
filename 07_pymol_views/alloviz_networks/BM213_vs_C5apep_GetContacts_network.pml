# BM213_vs_C5apep GetContacts node/edge delta
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
# Negative color: C5apep higher

# Key residue nodes. Cyan means positive state higher. Orange means comparator higher.
select all_key_nodes, none
select node_I91_58, receptor and resi 58
show sticks, node_I91_58
color delta_negative, node_I91_58
show spheres, node_I91_58 and name CA
set sphere_scale, 0.293, node_I91_58 and name CA
label node_I91_58 and name CA, "I91"
select all_key_nodes, all_key_nodes or node_I91_58
select node_W102_69, receptor and resi 69
show sticks, node_W102_69
color delta_positive, node_W102_69
show spheres, node_W102_69 and name CA
set sphere_scale, 0.564, node_W102_69 and name CA
label node_W102_69 and name CA, "W102"
select all_key_nodes, all_key_nodes or node_W102_69
select node_I116_83, receptor and resi 83
show sticks, node_I116_83
color delta_negative, node_I116_83
show spheres, node_I116_83 and name CA
set sphere_scale, 0.422, node_I116_83 and name CA
label node_I116_83 and name CA, "I116"
select all_key_nodes, all_key_nodes or node_I116_83
select node_M120_87, receptor and resi 87
show sticks, node_M120_87
color delta_positive, node_M120_87
show spheres, node_M120_87 and name CA
set sphere_scale, 0.498, node_M120_87 and name CA
label node_M120_87 and name CA, "M120"
select all_key_nodes, all_key_nodes or node_M120_87
select node_R134_101, receptor and resi 101
show sticks, node_R134_101
color delta_negative, node_R134_101
show spheres, node_R134_101 and name CA
set sphere_scale, 0.369, node_R134_101 and name CA
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
color delta_positive, node_Y222_189
show spheres, node_Y222_189 and name CA
set sphere_scale, 0.280, node_Y222_189 and name CA
label node_Y222_189 and name CA, "Y222"
select all_key_nodes, all_key_nodes or node_Y222_189
select node_F251_218, receptor and resi 218
show sticks, node_F251_218
color delta_negative, node_F251_218
show spheres, node_F251_218 and name CA
set sphere_scale, 0.265, node_F251_218 and name CA
label node_F251_218 and name CA, "F251"
select all_key_nodes, all_key_nodes or node_F251_218
select node_W255_222, receptor and resi 222
show sticks, node_W255_222
color delta_positive, node_W255_222
show spheres, node_W255_222 and name CA
set sphere_scale, 0.379, node_W255_222 and name CA
label node_W255_222 and name CA, "W255"
select all_key_nodes, all_key_nodes or node_W255_222
select node_Y258_225, receptor and resi 225
show sticks, node_Y258_225
color delta_positive, node_Y258_225
show spheres, node_Y258_225 and name CA
set sphere_scale, 0.345, node_Y258_225 and name CA
label node_Y258_225 and name CA, "Y258"
select all_key_nodes, all_key_nodes or node_Y258_225
select node_Y290_257, receptor and resi 257
show sticks, node_Y290_257
color delta_positive, node_Y290_257
show spheres, node_Y290_257 and name CA
set sphere_scale, 0.566, node_Y290_257 and name CA
label node_Y290_257 and name CA, "Y290"
select all_key_nodes, all_key_nodes or node_Y290_257
select node_N292_259, receptor and resi 259
show sticks, node_N292_259
color delta_neutral, node_N292_259
show spheres, node_N292_259 and name CA
set sphere_scale, 0.250, node_N292_259 and name CA
label node_N292_259 and name CA, "N292"
select all_key_nodes, all_key_nodes or node_N292_259
select node_N296_263, receptor and resi 263
show sticks, node_N296_263
color delta_negative, node_N296_263
show spheres, node_N296_263 and name CA
set sphere_scale, 0.496, node_N296_263 and name CA
label node_N296_263 and name CA, "N296"
select all_key_nodes, all_key_nodes or node_N296_263
select node_Y300_267, receptor and resi 267
show sticks, node_Y300_267
color delta_negative, node_Y300_267
show spheres, node_Y300_267 and name CA
set sphere_scale, 0.329, node_Y300_267 and name CA
label node_Y300_267 and name CA, "Y300"
select all_key_nodes, all_key_nodes or node_Y300_267

orient all_key_nodes
zoom all_key_nodes, 12

# Top delta edges shown as CA-CA dashed links.
select all_delta_edges, none
distance edge_1_R145_Y159, receptor and resi 145 and name CA, receptor and resi 159 and name CA
hide labels, edge_1_R145_Y159
color edge_negative, edge_1_R145_Y159
set dash_width, 4.65, edge_1_R145_Y159
select all_delta_edges, all_delta_edges or receptor and resi 145+159
distance edge_2_Y222_F251, receptor and resi 189 and name CA, receptor and resi 218 and name CA
hide labels, edge_2_Y222_F251
color edge_negative, edge_2_Y222_F251
set dash_width, 4.28, edge_2_Y222_F251
select all_delta_edges, all_delta_edges or receptor and resi 189+218
distance edge_3_I91_Q226, receptor and resi 91 and name CA, receptor and resi 226 and name CA
hide labels, edge_3_I91_Q226
color edge_negative, edge_3_I91_Q226
set dash_width, 3.94, edge_3_I91_Q226
select all_delta_edges, all_delta_edges or receptor and resi 91+226
distance edge_4_I91_F219, receptor and resi 91 and name CA, receptor and resi 219 and name CA
hide labels, edge_4_I91_F219
color edge_negative, edge_4_I91_F219
set dash_width, 3.92, edge_4_I91_F219
select all_delta_edges, all_delta_edges or receptor and resi 91+219
distance edge_5_F251_Q226, receptor and resi 218 and name CA, receptor and resi 226 and name CA
hide labels, edge_5_F251_Q226
color edge_negative, edge_5_F251_Q226
set dash_width, 3.87, edge_5_F251_Q226
select all_delta_edges, all_delta_edges or receptor and resi 218+226
distance edge_6_H68_K246, receptor and resi 68 and name CA, receptor and resi 246 and name CA
hide labels, edge_6_H68_K246
color edge_positive, edge_6_H68_K246
set dash_width, 3.86, edge_6_H68_K246
select all_delta_edges, all_delta_edges or receptor and resi 68+246
distance edge_7_I91_I116, receptor and resi 58 and name CA, receptor and resi 83 and name CA
hide labels, edge_7_I91_I116
color edge_negative, edge_7_I91_I116
set dash_width, 3.64, edge_7_I91_I116
select all_delta_edges, all_delta_edges or receptor and resi 58+83
distance edge_8_L94_V211, receptor and resi 94 and name CA, receptor and resi 211 and name CA
hide labels, edge_8_L94_V211
color edge_negative, edge_8_L94_V211
set dash_width, 3.63, edge_8_L94_V211
select all_delta_edges, all_delta_edges or receptor and resi 94+211
distance edge_9_L59_W102, receptor and resi 59 and name CA, receptor and resi 69 and name CA
hide labels, edge_9_L59_W102
color edge_positive, edge_9_L59_W102
set dash_width, 3.55, edge_9_L59_W102
select all_delta_edges, all_delta_edges or receptor and resi 59+69
distance edge_10_M120_Q226, receptor and resi 87 and name CA, receptor and resi 226 and name CA
hide labels, edge_10_M120_Q226
color edge_negative, edge_10_M120_Q226
set dash_width, 3.49, edge_10_M120_Q226
select all_delta_edges, all_delta_edges or receptor and resi 87+226
distance edge_11_Y159_R173, receptor and resi 159 and name CA, receptor and resi 173 and name CA
hide labels, edge_11_Y159_R173
color edge_negative, edge_11_Y159_R173
set dash_width, 3.48, edge_11_Y159_R173
select all_delta_edges, all_delta_edges or receptor and resi 159+173
distance edge_12_L59_Y290, receptor and resi 59 and name CA, receptor and resi 257 and name CA
hide labels, edge_12_L59_Y290
color edge_positive, edge_12_L59_Y290
set dash_width, 3.46, edge_12_L59_Y290
select all_delta_edges, all_delta_edges or receptor and resi 59+257
distance edge_13_R145_P150, receptor and resi 145 and name CA, receptor and resi 150 and name CA
hide labels, edge_13_R145_P150
color edge_positive, edge_13_R145_P150
set dash_width, 3.11, edge_13_R145_P150
select all_delta_edges, all_delta_edges or receptor and resi 145+150
distance edge_14_M120_Y258, receptor and resi 87 and name CA, receptor and resi 225 and name CA
hide labels, edge_14_M120_Y258
color edge_positive, edge_14_M120_Y258
set dash_width, 3.08, edge_14_M120_Y258
select all_delta_edges, all_delta_edges or receptor and resi 87+225
distance edge_15_Y88_Q226, receptor and resi 88 and name CA, receptor and resi 226 and name CA
hide labels, edge_15_Y88_Q226
color edge_negative, edge_15_Y88_Q226
set dash_width, 3.04, edge_15_Y88_Q226
select all_delta_edges, all_delta_edges or receptor and resi 88+226


set stick_radius, 0.16, all_key_nodes
show cartoon, receptor
zoom all_key_nodes, 10
# Optional rendering commands:
# ray 2400, 1800
# png BM213_vs_C5apep_GetContacts_network.png, dpi=300

