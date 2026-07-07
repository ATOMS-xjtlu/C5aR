#!/bin/bash

#SBATCH -J bm213_2        # Job name
#SBATCH --partition=gpu4090     # gpua800
#SBATCH --qos=4gpus
#SBATCH -N 1                    # Single node
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1            # 1 GPUs
#SBATCH -o %j.out         # Output result
#SBATCH -e %j.err          # Error output


ml load amber

parm complex-amber.top
trajin production.trj
trajin production1.trj
trajin production2.trj
#Hrmsdrmsd-dataset
rms rmsd-dataset first !@H=

#average-dataset
average crdset average-dataset

#1~30CA
#ca_matrix.gnuca_matrix:
matrix covar out ca_matrix.gnu name ca_matrix start 1 stop 20000 @CA

run
