#!/bin/bash
#SBATCH --job-name=yt
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=16
#SBATCH --time=10:00:00
#SBATCH --output=/home/fcastillo/yt.out
#SBATCH --mail-user=fabien.castillo@etu.unistra.fr
#SBATCH --mail-type=END,FAIL 

module purge
module load intelpython

mpirun -np 16 python /home/fcastillo/Cosmologicateur-Genial/Cosmologicateur.py

exit 0
