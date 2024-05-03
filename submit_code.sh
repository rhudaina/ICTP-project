#!/bin/bash -l
#SBATCH -A tra24_ictp_np
#SBATCH -p boost_usr_prod
#SBATCH --time 1:15:00       # format: HH:MM:SS
#SBATCH -N 1
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=1
#SBATCH --gpus-per-node=2
#SBATCH --mem-per-cpu=10000
#SBATCH --job-name=cell_project
#SBATCH --output=celljob_output.txt
#SBATCH --error=celljob_output.err

mkdir ./output

module load openmpi

echo -e "Running..."
mpicc -O3 main.c
mpirun -np 2 ./a.out 750 1000000

echo -e "Done!"
