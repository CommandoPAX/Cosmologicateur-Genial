#!/bin/bash
#SBATCH --job-name=yt
#SBATCH --nodes=1
#SBATCH --mem=100gb
#SBATCH --time=02:00:00
#SBATCH --output=/home/fcastillo/logs/yt.out

module purge
module load intelpython

python Cosmologicateur.py