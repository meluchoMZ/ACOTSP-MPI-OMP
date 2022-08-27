#!/bin/bash
#
#SBATCH --job-name=ft_aco
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -c 1

export PATH="$PATH:/home/miguel.blanco/mpi402/bin"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/home/miguel.blanco/mpi402/lib"
export I_MPI_FAULT_CONTINUE=1

# min pr2392 -> 378032
# min lin319

mpirun -np 16 ./acotsp -i pr2392.tsp -r 1 -t 1000 -x -o  378800
mpirun -np 16 ./acotsp -i pr2392.tsp -r 1 -t 1000 -x -o  378800
mpirun -np 16 ./acotsp -i pr2392.tsp -r 1 -t 1000 -x -o  378800
mpirun -np 16 ./acotsp -i pr2392.tsp -r 1 -t 1000 -x -o  378800
mpirun -np 16 ./acotsp -i pr2392.tsp -r 1 -t 1000 -x -o  378800
mpirun -np 16 ./acotsp -i pr2392.tsp -r 1 -t 1000 -x -o  378800
mpirun -np 16 ./acotsp -i pr2392.tsp -r 1 -t 1000 -x -o  378800
mpirun -np 16 ./acotsp -i pr2392.tsp -r 1 -t 1000 -x -o  378800
mpirun -np 16 ./acotsp -i pr2392.tsp -r 1 -t 1000 -x -o  378800
mpirun -np 16 ./acotsp -i pr2392.tsp -r 1 -t 1000 -x -o  378800
mpirun -np 16 ./acotsp -i pr2392.tsp -r 1 -t 1000 -x -o  378800
mpirun -np 16 ./acotsp -i pr2392.tsp -r 1 -t 1000 -x -o  378800
mpirun -np 16 ./acotsp -i pr2392.tsp -r 1 -t 1000 -x -o  378800
mpirun -np 16 ./acotsp -i pr2392.tsp -r 1 -t 1000 -x -o  378800
mpirun -np 16 ./acotsp -i pr2392.tsp -r 1 -t 1000 -x -o  378800
mpirun -np 16 ./acotsp -i pr2392.tsp -r 1 -t 1000 -x -o  378800
mpirun -np 16 ./acotsp -i pr2392.tsp -r 1 -t 1000 -x -o  378800
mpirun -np 16 ./acotsp -i pr2392.tsp -r 1 -t 1000 -x -o  378800
mpirun -np 16 ./acotsp -i pr2392.tsp -r 1 -t 1000 -x -o  378800
mpirun -np 16 ./acotsp -i pr2392.tsp -r 1 -t 1000 -x -o  378800
mpirun -np 16 ./acotsp -i pr2392.tsp -r 1 -t 1000 -x -o  378800
mpirun -np 16 ./acotsp -i pr2392.tsp -r 1 -t 1000 -x -o  378800
mpirun -np 16 ./acotsp -i pr2392.tsp -r 1 -t 1000 -x -o  378800
mpirun -np 16 ./acotsp -i pr2392.tsp -r 1 -t 1000 -x -o  378800
mpirun -np 16 ./acotsp -i pr2392.tsp -r 1 -t 1000 -x -o  378800
mpirun -np 16 ./acotsp -i pr2392.tsp -r 1 -t 1000 -x -o  378800
mpirun -np 16 ./acotsp -i pr2392.tsp -r 1 -t 1000 -x -o  378800
mpirun -np 16 ./acotsp -i pr2392.tsp -r 1 -t 1000 -x -o  378800
mpirun -np 16 ./acotsp -i pr2392.tsp -r 1 -t 1000 -x -o  378800
mpirun -np 16 ./acotsp -i pr2392.tsp -r 1 -t 1000 -x -o  378800
mpirun -np 16 ./acotsp -i pr2392.tsp -r 1 -t 1000 -x -o  378800
mpirun -np 16 ./acotsp -i pr2392.tsp -r 1 -t 1000 -x -o  378800
mpirun -np 16 ./acotsp -i pr2392.tsp -r 1 -t 1000 -x -o  378800
mpirun -np 16 ./acotsp -i pr2392.tsp -r 1 -t 1000 -x -o  378800
mpirun -np 16 ./acotsp -i pr2392.tsp -r 1 -t 1000 -x -o  378800
mpirun -np 16 ./acotsp -i pr2392.tsp -r 1 -t 1000 -x -o  378800
mpirun -np 16 ./acotsp -i pr2392.tsp -r 1 -t 1000 -x -o  378800
mpirun -np 16 ./acotsp -i pr2392.tsp -r 1 -t 1000 -x -o  378800
mpirun -np 16 ./acotsp -i pr2392.tsp -r 1 -t 1000 -x -o  378800
mpirun -np 16 ./acotsp -i pr2392.tsp -r 1 -t 1000 -x -o  378800
mpirun -np 16 ./acotsp -i pr2392.tsp -r 1 -t 1000 -x -o  378800
mpirun -np 16 ./acotsp -i pr2392.tsp -r 1 -t 1000 -x -o  378800
mpirun -np 16 ./acotsp -i pr2392.tsp -r 1 -t 1000 -x -o  378800
mpirun -np 16 ./acotsp -i pr2392.tsp -r 1 -t 1000 -x -o  378800
mpirun -np 16 ./acotsp -i pr2392.tsp -r 1 -t 1000 -x -o  378800
mpirun -np 16 ./acotsp -i pr2392.tsp -r 1 -t 1000 -x -o  378800
mpirun -np 16 ./acotsp -i pr2392.tsp -r 1 -t 1000 -x -o  378800
mpirun -np 16 ./acotsp -i pr2392.tsp -r 1 -t 1000 -x -o  378800
mpirun -np 16 ./acotsp -i pr2392.tsp -r 1 -t 1000 -x -o  378800
mpirun -np 16 ./acotsp -i pr2392.tsp -r 1 -t 1000 -x -o  378800
#mpiexec -n 16 --mca mpi_ft_enable=true ./acotsp -i pr2392.tsp -r 10 -t 100 -x -o 378032

