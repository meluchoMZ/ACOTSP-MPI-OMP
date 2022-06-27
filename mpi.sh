#!/bin/bash
#
#SBATCH --job-name=ft_aco
#SBATCH -N 1
#SBATCH -n 4 
#SBATCH -c 1

export PATH="$PATH:/home/miguel.blanco/mpi402/bin"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/home/miguel.blanco/mpi402/lib"
export I_MPI_FAULT_CONTINUE=1
mpirun -np 4 ./acotsp -i pr2392.tsp -r 3 -t 10 -x -o 378032
#mpiexec -n 16 --mca mpi_ft_enable=true ./acotsp -i pr2392.tsp -r 10 -t 100 -x -o 378032

