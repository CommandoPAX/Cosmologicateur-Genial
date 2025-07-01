import os
import sys
from astropy.io import fits
import numpy as np
import time

pre = "/data100/fcastillo/RESULT/"
snapshots = ["benchM","NG_F500","G_m500","NG_F500_m500","NG_Fminus500","NG_Fminus500_m500","G_ViVi", "NG_ViVi","NG_Fminus500_ViVi","NEDE","NsPNG_EDE_F500","NsPNG_EDE_F1833", "NBM"]

indices_hdm = [0,1,4,6,7,8,9,10,11]
indices_z = [6,9]

#indices_hdm = [0,2,6]
#indices_hdm = [0,1,4,6,9]

for n in [12]: # 7
    if n <= 8 : redshifts = [2,4]
    else : redshifts = indices_z
    for i in redshifts:
        #input_ = pre + snapshots[n]+"/"+str(i)+"_densite"


        fichier = open(f"/home/fcastillo/bash/extrema_{n}_{i}.sh","w")
        fichier.write(f"""#!/bin/bash
#SBATCH --job-name=extrema_{n}_{i}
#SBATCH --nodes=1
#SBATCH --time=48:00:00
#SBATCH --output=/home/fcastillo/logs/extrema_{n}_{i}.out
#SBATCH --partition=pscomp
#SBATCH --ntasks-per-node=128
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
python /home/fcastillo/Cosmologicateur-Genial/extremateur.py {n} {i}
""")
            
        fichier.write("exit 0")
        fichier.close()
        #os.system(f"sbatch extrema_{n}_{i}.sh")
