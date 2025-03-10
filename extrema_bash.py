import os
import sys
from astropy.io import fits
import numpy as np
import time

pre = "/data100/fcastillo/RESULT/"
snapshots = ["benchM","NG_F500","G_m500","NG_F500_m500","NG_Fminus500","NG_Fminus500_m500"]

for n in range(6):
    for i in range(5):

        input_ = pre + snapshots[n]+"/"+str(i)+"_densite"


        fichier = open(f"extrema_{n}_{i}_2.sh","w")
        fichier.write(f"""#!/bin/bash
#SBATCH --job-name=extrema_{n}_{i}
#SBATCH --nodes=1
#SBATCH --mem=500gb
#SBATCH --time=02:00:00
#SBATCH --output=/home/fcastillo/logs/extrema_{n}_{i}.out
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


        fichier.write(f"""
python extremateur.py {n} {i}
""")
            
        fichier.write("exit 0")
        fichier.close()
        os.system(f"sbatch extrema_{n}_{i}.sh")
