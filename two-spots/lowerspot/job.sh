#!/bin/bash
#SBATCH -N 1 -n 64
#SBATCH --partition=neutron-star
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --job-name=pix300lo
#SBATCH --mail-type=ALL
#SBATCH --array=0-100
#SBATCH --time=00-01:30:00
#SBATCH -o outfiles/out_%a # STDOUT 
#SBATCH -e errfiles/err_%a # STDERR

OMP_STACKSIZE=20m
export OMP_STACKSIZE

./RAPTOR model.in /zfs/helios/filer0/pdas1/ray_tracing/RAPTOR/data_files/inc60/data0508.dat 1 $SLURM_ARRAY_TASK_ID
