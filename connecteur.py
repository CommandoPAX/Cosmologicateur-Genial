import os
import sys
from astropy.io import fits
import numpy as np
import time

pre = "/data100/fcastillo/RESULT/"

snapshots = ["NEDE","NsPNG_F500","NsPNG_F1000","NsPNG_F1833","NsPNG_EDE_F500","NsPNG_EDE_F1000","NsPNG_EDE_F1833"]

z= [15,12, 10, 8, 5,3,1,0.5,0.25,0]
indices_z = [5,6,8,9]


for n in range(7):
    for i in indices_z:

        #input_ = pre + snapshots[n]+"/"+str(i)+"_densite"


        fichier = open(f"connect_{n}_{i}.sh","w")
        fichier.write(f"""#!/bin/bash
#SBATCH --job-name=fil_{n}_{i}
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
python connect_filaments.py {n} {i}
""")
            
        fichier.write("exit 0")
        fichier.close()
        #os.system(f"sbatch connect_{n}_{i}.sh")
