import os
import sys
from astropy.io import fits
import numpy as np
import time

pre = "/data100/fcastillo/RESULT/"
snapshots = ["benchM","NG_F500","G_m500","NG_F500_m500","NG_Fminus500","NG_Fminus500_m500", "G_ViVi", "NG_ViVi", "NG_Fminus500_ViVi"]

for n in range(6,9):
    for i in range(1,5):

        input_ = pre + snapshots[n]+"/"+str(i)+"_densite"


        fichier = open(f"connect_{n}_{i}.sh","w")
        fichier.write(f"""#!/bin/bash
#SBATCH --job-name=halo_{n}_{i}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100gb
#SBATCH --time=2:00:00
#SBATCH --output=/home/fcastillo/connect_{n}_{i}.out
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


        fichier.write(f"""
python connecteur_halos.py {n} {i}
""")
            
        fichier.write("exit 0")
        fichier.close()
        os.system(f"sbatch connect_{n}_{i}.sh")
