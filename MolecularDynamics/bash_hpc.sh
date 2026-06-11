#!/bin/bash -l
#
#SBATCH --gres=gpu:a40:1
#SBATCH --time=01:00:00
#SBATCH --export=NONE

unset SLURM_EXPORT_ENV
module load nvhpc

rm -r -f out_cuda/
mkdir -p out_cuda/

make cuda
rm -f particle_cuda.tar.gz
tar -cvf particle_cuda.tar.gz *.vtu *.pvd
rm -r *.vtu *.pvd 
rm -r -f out_cuda/