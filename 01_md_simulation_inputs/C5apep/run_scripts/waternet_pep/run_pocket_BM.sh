#!/bin/bash

#SBATCH -J pep_pocket
#SBATCH --partition=cpu8358
#SBATCH --qos=52cores
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH -o %j.out
#SBATCH -e %j.err

echo "==== Pocket water job started: `date` ===="

# 
WATERNET_ENV="/gpfs/work/pha/linhaozheng22/waternet"
PYTHON="$WATERNET_ENV/bin/python"

echo "Using python: $PYTHON"

# topology
TOP="../md-initial_WAT.pdb"

# xtc 
TRJ1="../4000frame_WAT_pep1.xtc"
TRJ2="../second/4000frame_WAT_pep2.xtc"
TRJ3="../third/4000frame_WAT_pep3.xtc"

OUTPREFIX="pep_pocket"

$PYTHON pocket_water_BM.py \
    -t $TOP \
    -x $TRJ1 $TRJ2 $TRJ3 \
    --resY6 191 \
    --resY7 269 \
    --radius 4.0 \
    --stride 1 \
    -o $OUTPREFIX

echo "==== Pocket water job finished: `date` ===="

