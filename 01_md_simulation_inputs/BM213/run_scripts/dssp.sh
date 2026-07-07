#!/bin/bash

#SBATCH -J dssp         # Job name
#SBATCH --partition=gpu4090     # gpua800
#SBATCH --qos=4gpus
#SBATCH -N 1                    # Single node
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1            # 1 GPUs
#SBATCH -o %j.out         # Output result
#SBATCH -e %j.err          # Error output

ml load gromacs/2023.2-gcc-9.5.0-jzxesel


gmx_mpi do_dssp -s md-initial.pdb -f md-nW.xtc -n dssp.ndx -o dssp.xpm
