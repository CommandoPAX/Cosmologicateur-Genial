import os
import sys
from astropy.io import fits
import numpy as np
import time

pre = "/data100/fcastillo/RESULT/"
snapshots = ["benchM","NG_F500","G_m500","NG_F500_m500","NG_Fminus500","NG_Fminus500_m500"]

fichier = open("mse.sh","w")
fichier.write(f"""#!/bin/bash
#SBATCH --job-name=disperse
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --ntasks-per-node=32
#SBATCH --mem=500gb
#SBATCH --time=24:00:00
#SBATCH --output=/home/fcastillo/disperse.out
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
module load disperse/0.9.24
""")

for n in range(6):
    for i in range(5):

        input_ = pre + snapshots[n]+"/"+str(i)+"_densite"

        fichier.write(f"""

        /softs/disperse/0.9.24/bin/mse {input_}_0.fits -nthreads 32 -cut 1 -upSkl -manifolds -outName {pre+snapshots[n]+"/"+str(i)+"_densite_0.fits"} -periodicity 0 -forceLoops
        /softs/disperse/0.9.24/bin/skelconv {pre+snapshots[n]+"/"+str(i)+"_densite_0.fits_c1.up.NDskl"} -smooth 1 -outName {pre+snapshots[n]+"/"+str(i)+"_densite_0.fits_c1.up.NDskl"} -to NDskl_ascii

        exit 0""")
            
        fichier.close()
os.system("sbatch mse.sh")
