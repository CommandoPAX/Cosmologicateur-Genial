#!/bin/bash
#SBATCH --job-name=pow
#SBATCH --partition=comp,pscomp
#SBATCH --nodes=1
#SBATCH --mem=10gb
#SBATCH --time=48:00:00
#SBATCH --output=/home/fcastillo/logs/plot.out

module purge
module load intelpython

python /home/fcastillo/Cosmologicateur-Genial/Cosmologicateur.py
python /home/fcastillo/Cosmologicateur-Genial/Analysateur_genial_simplifie.py
