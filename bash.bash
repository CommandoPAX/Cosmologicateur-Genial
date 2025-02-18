#!/bin/bash
#SBATCH --job-name=yt
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --ntasks-per-node=32
#SBATCH --time=10:00:00
#SBATCH --output=/home/fcastillo/yt.out
#SBATCH --mail-user=fabien.castillo@etu.unistra.fr
#SBATCH --mail-type=END,FAIL 

module purge
module load intelpython
module load inteloneapi/2025.0.1

mpirun -np 32 python /home/fcastillo/Cosmologicateur-Genial/Cosmologicateur.py

exit 0
