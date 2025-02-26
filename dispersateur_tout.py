import os
import sys
from astropy.io import fits
import numpy as np
import time

pre = "/data100/fcastillo/RESULT/"
snapshots = ["benchM","NG_F500","G_m500","NG_F500_m500","NG_Fminus500","NG_Fminus500_m500"]


for n in range(6):
    for i in range(5):
        fichier = open("mse_"+str(n)+"_"+str(i)+".sh","w")

        input_ = pre + snapshots[n]+"/"+str(i)+"_densite"

        fichier.write(f"""#!/bin/bash
#SBATCH --job-name=disperse_{n}_{i}
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --ntasks-per-node=32
#SBATCH --mem=200gb
#SBATCH --time=2:00:00
#SBATCH --output=/home/fcastillo/disperse_{n}_{i}.out
#SBATCH --mail-user=fabien.castillo@etu.unistra.fr
#SBATCH --mail-type=ALL
#SBATCH --partition=pscomp""")
        fichier.write("""

module () {
  eval $(/usr/bin/modulecmd bash $*)
}

source /etc/profile.d/modules.sh
""")

        fichier.write(f"""
module purge
module load disperse/0.9.24

/softs/disperse/0.9.24/bin/mse {input_}_0.fits -nthreads 32 -upSkl -manifolds -outName {pre+snapshots[n]+"/"+str(i)+"_densite_0_time.fits"} -periodicity 0 -forceLoops
/softs/disperse/0.9.24/bin/skelconv {pre+snapshots[n]+"/"+str(i)+"_densite_0_time.fits_c1.up.NDskl"} -smooth 1 -outName {pre+snapshots[n]+"/"+str(i)+"_densite_0_time.fits_c1.up.NDskl"} -to NDskl_ascii

exit 0""")
        
        fichier.close()
        os.system("sbatch mse_"+str(n)+"_"+str(i)+".sh")
        time.sleep(20*60)