#!/bin/bash
#SBATCH --job-name=lisse
#SBATCH --nodes=1
#SBATCH --mem=100gb
#SBATCH --time=02:00:00
#SBATCH --output=/home/fcastillo/logs/lisse.out

module purge
module load intelpython

python dispersateur_tout.py