#!/bin/bash

#SBATCH -J c5apep_second        # Job name
#SBATCH --partition=gpu4090     # gpua800
#SBATCH --qos=4gpus
#SBATCH -N 1                    # Single node
#SBATCH --ntasks-per-node=1     
#SBATCH --cpus-per-task=4       
#SBATCH --gres=gpu:1            # 1 GPUs
#SBATCH -o %j.out         # Output result
#SBATCH -e %j.err          # Error output


ml load amber

#istart=1
#iend=3

#I=$istart
#while [ $I -le $iend ]
#do
#  suffix=`printf %03d $I`
#  suffix2=`printf %03d $((I-1))`
#  pmemd.cuda -O -i prod.in -o prod$suffix.out -p wbox.prmtop -c prod$suffix2.rst -x prod$suffix.nc -r prod$suffix.rst
#  I=$((I+1))
#done
pmemd.cuda -O -i ../7_production.mdin -p ../input.parm7 -c ../6_equilibration.rst -o production.out -x production.trj -r production.rst -ref ../6_equilibration.rst
