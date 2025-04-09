#!/bin/bash
#SBATCH --job-name=plot
#SBATCH --partition=comp,pscomp
#SBATCH --nodes=1
#SBATCH --mem=10gb
#SBATCH --time=00:01:00
#SBATCH --output=/home/fcastillo/logs/plot.out

module purge
module load intelpython

python plot_minkowski2.py
