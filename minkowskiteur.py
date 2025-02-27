import os
import sys
from astropy.io import fits
import numpy as np
import time

pre = "/data100/fcastillo/RESULT/"
snapshots = ["benchM","NG_F500","G_m500","NG_F500_m500","NG_Fminus500","NG_Fminus500_m500"]

fichier = open("minkowski.sh","w")
fichier.write(f"""#!/bin/bash
#SBATCH --job-name=minkowski
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --ntasks-per-node=32
#SBATCH --mem=500gb
#SBATCH --time=48:00:00
#SBATCH --output=/home/fcastillo/minkowski.out
#SBATCH --mail-user=fabien.castillo@etu.unistra.fr
#SBATCH --mail-type=ALL
#SBATCH --partition=pscomp
""")
fichier.write("""

module () {
eval $(/usr/bin/modulecmd bash $*)
}

source /etc/profile.d/modules.sh
""")

fichier.write(f"""
module purge
module load intelpython
module load inteloneapi/2025.0.1
""")

for n in range(6):
    for i in range(5):

        input_ = pre + snapshots[n]+"/"+str(i)+"_densite"

        fichier.write(f"""
python minkowski.py {n} {i}
""")
            
fichier.write("exit 0")
fichier.close()
os.system("sbatch minkowski.sh")
