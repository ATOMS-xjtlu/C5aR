# β-GPR1  C5AR 

## 

- `article/analysis/sommd/pymol/formal/β-GPR1.pdb`
- chain `R` = GPR1/GPCRchain `A` = β-arrestin 
- C5AR representative`article/analysis/sommd/pymol/formal/cluster_01_BM213_rep3.pdb`  `article/analysis/sommd/pymol/formal/cluster_02_C5apep_rep3.pdb`
- getcontacts`article/analysis/sommd/pymol/formal/interface_contacts/beta_gpr1_getcontacts_chainR_chainA.tsv` `--sele 'chain R' --sele2 'chain A' --itypes all --distout`
- GPR1→C5AR  chain R  C5AR chain X  BLOSUM62  C5AR docking  C5AR  β-arrestin

## 

- getcontacts  β-GPR1  GPR1  `18`  active/core  `10` 
- β-arrestin  `13`  active/core  `10` 
- 4.5 Å heavy-atom / getcontacts  active residues
- C5AR docking PDB  chain `X`  1-281C5AR / =  + `33` offset  activation-core mapping  101  C5AR R134

## HADDOCK  active/passive 

HADDOCK  active residues  AIR HDOCK  C5AR  C5AR  docking PDB  HADDOCK 

### C5AR receptor

| role | C5AR  | docking PDB  |  | evidence |
|---|---|---|---|---:|
| active | P141 | X:PRO108 | GPR1 R:PRO142 | vdw:13; n=13; min=3.113 A |
| active | I142 | X:ILE109 | GPR1 R:VAL143 | vdw:3; n=3; min=2.974 A |
| active | T240 | X:THR207 | GPR1 R:HIS244 | vdw:8; n=8; min=2.977 A |
| active | V243 | X:VAL210 | GPR1 R:THR247 | vdw:2; n=2; min=3.121 A |
| active | G304 | X:GLY271 | GPR1 R:SER308 | vdw:3; n=3; min=3.277 A |
| active | Q305 | X:GLN272 | GPR1 R:LYS309 | sb:2;vdw:3; n=5; min=3.194 A |
| passive | T69 | X:THR36 | GPR1 R:THR70 | vdw:3; n=3; min=2.776 A |
| passive | R134 | X:ARG101 | GPR1 R:HIS135 | vdw:4; n=4; min=3.554 A |
| passive | F139 | X:PHE106 | GPR1 R:ILE140 | vdw:1; n=1; min=3.250 A |
| passive | A234 | X:ALA201 | GPR1 R:ILE238 | vdw:6; n=6; min=3.582 A |

### β-arrestin ligand

| role | β-arrestin  | HADDOCK  | evidence |
|---|---|---|---:|
| active | A:ASP67 | A:ASP67 | sb:2;vdw:2; n=4; min=3.214 A |
| active | A:LEU71 | A:LEU71 | vdw:11; n=11; min=2.977 A |
| active | A:LEU73 | A:LEU73 | vdw:6; n=6; min=3.412 A |
| active | A:PHE244 | A:PHE244 | vdw:7; n=7; min=3.582 A |
| active | A:TYR249 | A:TYR249 | vdw:13; n=13; min=2.974 A |
| passive | A:GLU66 | A:GLU66 | vdw:4; n=4; min=2.776 A |
| passive | A:ASP69 | A:ASP69 | vdw:3; n=3; min=3.277 A |
| passive | A:VAL70 | A:VAL70 | vdw:3; n=3; min=3.511 A |
| passive | A:ILE241 | A:ILE241 | vdw:3; n=3; min=3.757 A |
| passive | A:ASN245 | A:ASN245 | vdw:1; n=1; min=3.250 A |

###  HADDOCK residue lists

- C5AR active`P141, I142, T240, V243, G304, Q305`
- C5AR active`X:PRO108, X:ILE109, X:THR207, X:VAL210, X:GLY271, X:GLN272`
- C5AR passive`T69, R134, F139, A234`
- C5AR passive`X:THR36, X:ARG101, X:PHE106, X:ALA201`
- β-arrestin active`A:ASP67, A:LEU71, A:LEU73, A:PHE244, A:TYR249`
- β-arrestin passive`A:GLU66, A:ASP69, A:VAL70, A:ILE241, A:ASN245`

## HDOCK/

### C5AR receptor site residues

- `cluster_01_BM213_rep3` core/recommended`X:36, X:101, X:106, X:108, X:109, X:201, X:207, X:210, X:271, X:272`
- `cluster_01_BM213_rep3` all getcontacts-mapped`X:34, X:36, X:42, X:101, X:104, X:106, X:107, X:108, X:109, X:117, X:189, X:201, X:202, X:207, X:210, X:267, X:271, X:272`
- `cluster_01_BM213_rep3` 4.5 Å /passive`X:33, X:34, X:36, X:38, X:39, X:42, X:101, X:104, X:106, X:107, X:108, X:109, X:115, X:116, X:117, X:189, X:200, X:201, X:202, X:206, X:207, X:210, X:211, X:267, X:270, X:271, X:272`
- `cluster_02_C5apep_rep3` core/recommended`X:36, X:101, X:106, X:108, X:109, X:201, X:207, X:210, X:271, X:272`
- `cluster_02_C5apep_rep3` all getcontacts-mapped`X:34, X:36, X:42, X:101, X:104, X:106, X:107, X:108, X:109, X:117, X:189, X:201, X:202, X:207, X:210, X:267, X:271, X:272`
- `cluster_02_C5apep_rep3` 4.5 Å /passive`X:33, X:34, X:36, X:38, X:39, X:42, X:101, X:104, X:106, X:107, X:108, X:109, X:115, X:116, X:117, X:189, X:200, X:201, X:202, X:206, X:207, X:210, X:211, X:267, X:270, X:271, X:272`

### β-arrestin ligand site residues

- core/recommended`A:66, A:67, A:69, A:70, A:71, A:73, A:241, A:244, A:245, A:249`
- all getcontacts`A:63, A:66, A:67, A:69, A:70, A:71, A:73, A:134, A:136, A:241, A:244, A:245, A:249`
- 4.5 Å /passive`A:63, A:65, A:66, A:67, A:68, A:69, A:70, A:71, A:72, A:73, A:77, A:134, A:136, A:241, A:244, A:245, A:249`

## GPR1  C5AR 

- `cluster_01_BM213_rep3`C5AR residues=281aligned columns=273sequence identity=0.326score=405.5
- `cluster_02_C5apep_rep3`C5AR residues=281aligned columns=273sequence identity=0.326score=405.5

| GPR1 getcontacts residue | evidence | C5AR/BM213 mapped | C5AR/C5apep mapped | use |
|---|---:|---|---|---|
| R:LYS68 | vdw:1; n=1; min=3.897 Å | X:LYS34 | X:LYS34 | broader/passive |
| R:THR70 | vdw:3; n=3; min=2.776 Å | X:THR36 | X:THR36 | active/core |
| R:PHE76 | vdw:1; n=1; min=3.511 Å | X:PHE42 | X:PHE42 | broader/passive |
| R:HIS135 | vdw:4; n=4; min=3.554 Å | X:ARG101 | X:ARG101 | active/core |
| R:HIS138 | vdw:1; n=1; min=3.642 Å | X:LEU104 | X:LEU104 | broader/passive |
| R:ILE140 | vdw:1; n=1; min=3.250 Å | X:PHE106 | X:PHE106 | active/core |
| R:HIS141 | vdw:2; n=2; min=3.503 Å | X:LYS107 | X:LYS107 | broader/passive |
| R:PRO142 | vdw:13; n=13; min=3.113 Å | X:PRO108 | X:PRO108 | active/core |
| R:VAL143 | vdw:3; n=3; min=2.974 Å | X:ILE109 | X:ILE109 | active/core |
| R:LEU151 | vdw:1; n=1; min=3.831 Å | X:ALA117 | X:ALA117 | broader/passive |
| R:TYR226 | vdw:2; n=2; min=3.573 Å | X:TYR189 | X:TYR189 | broader/passive |
| R:ILE238 | vdw:6; n=6; min=3.582 Å | X:ALA201 | X:ALA201 | active/core |
| R:LEU239 | vdw:1; n=1; min=3.670 Å | X:THR202 | X:THR202 | broader/passive |
| R:HIS244 | vdw:8; n=8; min=2.977 Å | X:THR207 | X:THR207 | active/core |
| R:THR247 | vdw:2; n=2; min=3.121 Å | X:VAL210 | X:VAL210 | active/core |
| R:TYR304 | vdw:2; n=2; min=3.534 Å | X:TYR267 | X:TYR267 | broader/passive |
| R:SER308 | vdw:3; n=3; min=3.277 Å | X:GLY271 | X:GLY271 | active/core |
| R:LYS309 | sb:2;vdw:3; n=5; min=3.194 Å | X:GLN272 | X:GLN272 | active/core |

## β-arrestin 

| β-arrestin residue | evidence | partners | use |
|---|---:|---|---|
| A:TYR63 | vdw:2; n=2; min=3.503 Å | R:HIS141 | broader/passive |
| A:GLU66 | vdw:4; n=4; min=2.776 Å | R:THR70;R:LYS309 | active/core |
| A:ASP67 | sb:2;vdw:2; n=4; min=3.214 Å | R:LYS309 | active/core |
| A:ASP69 | vdw:3; n=3; min=3.277 Å | R:SER308 | active/core |
| A:VAL70 | vdw:3; n=3; min=3.511 Å | R:PHE76;R:HIS135;R:HIS138 | active/core |
| A:LEU71 | vdw:11; n=11; min=2.977 Å | R:HIS135;R:TYR226;R:HIS244;R:THR247;R:TYR304 | active/core |
| A:LEU73 | vdw:6; n=6; min=3.412 Å | R:HIS135;R:HIS244 | active/core |
| A:GLU134 | vdw:1; n=1; min=3.831 Å | R:LEU151 | broader/passive |
| A:THR136 | vdw:1; n=1; min=3.897 Å | R:LYS68 | broader/passive |
| A:ILE241 | vdw:3; n=3; min=3.757 Å | R:PRO142 | active/core |
| A:PHE244 | vdw:7; n=7; min=3.582 Å | R:ILE238;R:LEU239 | active/core |
| A:ASN245 | vdw:1; n=1; min=3.250 Å | R:ILE140 | active/core |
| A:TYR249 | vdw:13; n=13; min=2.974 Å | R:PRO142;R:VAL143 | active/core |

## 

| GPR1 | β-arrestin | evidence |
|---|---|---:|
| R:PRO142 | A:TYR249 | vdw:10; n=10; min=3.113 Å |
| R:ILE238 | A:PHE244 | vdw:6; n=6; min=3.582 Å |
| R:HIS244 | A:LEU73 | vdw:5; n=5; min=3.412 Å |
| R:LYS309 | A:ASP67 | sb:2;vdw:2; n=4; min=3.214 Å |
| R:THR70 | A:GLU66 | vdw:3; n=3; min=2.776 Å |
| R:VAL143 | A:TYR249 | vdw:3; n=3; min=2.974 Å |
| R:HIS244 | A:LEU71 | vdw:3; n=3; min=2.977 Å |
| R:SER308 | A:ASP69 | vdw:3; n=3; min=3.277 Å |
| R:PRO142 | A:ILE241 | vdw:3; n=3; min=3.757 Å |
| R:THR247 | A:LEU71 | vdw:2; n=2; min=3.121 Å |
| R:HIS141 | A:TYR63 | vdw:2; n=2; min=3.503 Å |
| R:TYR304 | A:LEU71 | vdw:2; n=2; min=3.534 Å |

##  docking 

- β-arrestin `article/analysis/sommd/pymol/formal/beta_arrestin_chainA_from_beta_GPR1.pdb`
- residue-level `article/analysis/sommd/pymol/formal/interface_contacts/beta_gpr1_getcontacts_interface_residues.csv`
- contact-pair `article/analysis/sommd/pymol/formal/interface_contacts/beta_gpr1_getcontacts_interface_pairs.csv`
- GPR1→C5AR `article/analysis/sommd/pymol/formal/interface_contacts/gpr1_interface_to_c5ar_mapping.csv`
- HADDOCK active/passive `article/analysis/sommd/pymol/formal/interface_contacts/haddock_active_passive_numbering.csv`
- BM213 HADDOCK receptor active `article/analysis/sommd/pymol/formal/interface_contacts/haddock_active_rsite_cluster_01_BM213_rep3.txt`
- C5apep HADDOCK receptor active `article/analysis/sommd/pymol/formal/interface_contacts/haddock_active_rsite_cluster_02_C5apep_rep3.txt`
- β-arrestin HADDOCK ligand active `article/analysis/sommd/pymol/formal/interface_contacts/haddock_active_lsite_beta_arrestin.txt`
- BM213 HDOCK receptor core `article/analysis/sommd/pymol/formal/interface_contacts/hdock_rsite_cluster_01_BM213_rep3_core.txt`
- C5apep HDOCK receptor core `article/analysis/sommd/pymol/formal/interface_contacts/hdock_rsite_cluster_02_C5apep_rep3_core.txt`
- β-arrestin HDOCK ligand core `article/analysis/sommd/pymol/formal/interface_contacts/hdock_lsite_beta_arrestin_core.txt`

##  HDOCK/HADDOCK 

- HDOCK C5AR receptor  core/recommended `rsite`β-arrestin ligand  core/recommended `lsite` top poses  all getcontacts 
- HADDOCK  active residues passive residues getcontacts  AIR
- BM213  C5apep  representative  C5AR  docking  β-arrestin  C5AR  intracellular pocket 
-  GPR1-β-arrestin  β-arrestin-compatible docking  C5AR  β-arrestin 
