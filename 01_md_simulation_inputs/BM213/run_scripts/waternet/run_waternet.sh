#!/bin/bash

#SBATCH -J bm213_waternet         # Job name
#SBATCH --partition=cpu8358     # gpua800
#SBATCH --qos=52cores
#SBATCH -N 1                    # Single node
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH -o %j.out         # Output result
#SBATCH -e %j.err          # Error output

echo "===== Job started at: `date` ====="

# 
WATERNET_ENV="/gpfs/work/pha/linhaozheng22/waternet"

PYTHON="$WATERNET_ENV/bin/python"

echo "Using python: $PYTHON"

# 
TOP="../md-initial_WAT.pdb"
TRJ1="../4500frame_WAT_bm1.xtc"
TRJ2="../second/4000frame_WAT_bm2.xtc"
TRJ3="../third/4000frame_WAT_bm3.xtc"

OUTPREFIX="bm213_water"

# 
$PYTHON water_bridges_pep.py \
    -t $TOP \
    -x $TRJ1 $TRJ2 $TRJ3 \
    --cutoff 4 \
    --stride 1 \
    -o $OUTPREFIX

echo "===== Job finished at: `date` ====="

