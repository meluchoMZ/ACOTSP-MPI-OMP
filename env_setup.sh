#!/bin/bash
## Environment setup
PARAMS=1

if [ $# -eq $PARAMS ]
then
	if [ $1 == "--start" ]
	then
		echo "$(date): loading kernel modules"
		echo "module load intel_parallel_studio_xe"
		module load intel_parallel_studio_xe
		echo "module load mkl"
		module load mkl
		echo "module load mpi"
		module load mpi
		module list
	else 
		if [ $1 == "--end" ]
		then
			echo "$(date): unloading kernel modules"
			echo "module unload mpi"
			module unload mpi
			echo "module unload mkl"
			module unload mkl
			echo "module unload intel_parallel_studio_xe"
			module unload intel_parallel_studio_xe
			echo "Loaded modules:"
			module list
		else
			echo "Error: unknown option '$1'"
			echo "Usage:"
			echo "	--start -> loads kernel modules"
			echo "	--end -> unloads kernel modules"
		fi
	fi
else
	echo "Error: syntax error"
	echo "Usage:"
	echo "	-start -> loads kernel modules"
	echo "	-end -> unloads kernel modules"
fi

