# C5aR1 Candidate Activation-Route Views

This folder contains direct-view structure visualizations for candidate BM213 and C5apep communication routes.

Use the HTML files for quick browser inspection. They embed the representative receptor PDB but load 3Dmol.js from `https://3Dmol.csb.pitt.edu`.

Use the PML files in Windows PyMOL for manuscript-grade ray-traced images.

| system | HTML | PyMOL |
|---|---|---|
| BM213 | `BM213_activation_route.html` | `BM213_activation_route.pml` |
| C5apep | `C5apep_activation_route.html` | `C5apep_activation_route.pml` |

Legend files: `activation_route_legend.png` and `activation_route_legend.svg`.

Interpretation boundary: these are candidate communication routes inferred from AlloViz-derived residue networks and NetworkX shortest paths. They should be combined with microswitch/FES/structural evidence before being described as activation mechanisms.

Route definitions are recorded in `activation_routes_summary.csv` and `activation_route_manifest.json`.
