#!/bin/bash
#SBATCH --job-name=plot
#SBATCH --nodes=1
#SBATCH --mem=100gb
#SBATCH --time=00:02:00
#SBATCH --output=/home/fcastillo/logs/plot.out

module purge
module load intelpython

python plot_threshold_2.py