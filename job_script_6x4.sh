#!/bin/bash
#
#SBATCH --job-name=test
#SBATCH -N 2 # (nodes)
#SBATCH -n 16 # (processes)

export HOME=/home/miguel.blanco/git/ACOTSP-MPI-OMP
export PATH="$PATH:/home/miguel.blanco/mpi5/bin"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH/home/miguel.blanco/mpi5/lib"

cd $HOME

mpirun -np $SLURM_NTASKS ./acotsp -i pr2392.tsp -r 10 -t 100 -x -o 378032
