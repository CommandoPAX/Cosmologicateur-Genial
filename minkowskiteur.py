import os
import sys
from astropy.io import fits
import numpy as np
import time

#pre = "/data100/fcastillo/RESULT/"
#snapshots = ["benchM","NG_F500","G_m500","NG_F500_m500","NG_Fminus500","NG_Fminus500_m500"]
pre = "/data77/stahl/Scale/Nb/WDM/ViVi/"
snapshots = ["G_ViVi","NG_Fminus500_ViVi","NG_ViVi"]
snapshots=["NBM"]

for n in range(1):
    if True :#for i in range(5):


        fichier = open(f"../bash/minkowski_{n}.sh","w")
        fichier.write(f"""#!/bin/bash
#SBATCH --job-name=minkowski_{n}
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --ntasks-per-node=128
#SBATCH --mem=100gb
#SBATCH --time=48:00:00
#SBATCH --output=/home/fcastillo/logs/minkowski_{n}.out
#SBATCH --partition=pscomp
#SBATCH --nodelist=i31,i32
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


        fichier.write(f"""
python /home/fcastillo/Cosmologicateur-Genial/minkowski.py {n}
""")
            
        fichier.write("exit 0")
        fichier.close()
        #os.system(f"sbatch ../bash/minkowski_{n}_{i}.sh")
