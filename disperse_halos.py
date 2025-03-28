import os
import sys
from astropy.io import fits
import numpy as np
import time

pre = "/data100/fcastillo/RESULT/"
snapshots = ["benchM","NG_F500","G_m500","NG_F500_m500","NG_Fminus500","NG_Fminus500_m500", "G_ViVi", "NG_ViVi", "NG_Fminus500_ViVi"]
for n in range(9):
    fichier = open(f"../bash/mse_{n}.sh","w")
    fichier.write(f"""#!/bin/bash
#SBATCH --job-name=skl2_{n}
#SBATCH --mem=100gb
#SBATCH --time=2:00:00
#SBATCH --output=/home/fcastillo/logs/skl2_{n}.out
""")
    fichier.write("""
module () {
  eval $(/usr/bin/modulecmd bash $*)
}

source /etc/profile.d/modules.sh
module purge
module load disperse/0.9.24""")

    for i in range(1,5):

        fichier.write(f"""
path_dat="{pre+snapshots[n]}/"
path="{pre+snapshots[n]}/"

file="{i}_halos"

persistence=1
smoothing=1

/softs/disperse/0.9.24/bin/delaunay_3D $path_dat$file".fits" -outName $path$file 

/softs/disperse/0.9.24/bin/mse $path_dat$file".NDnet" -nsig $persistence -upSkl -outName $path$file -forceLoops
/softs/disperse/0.9.24/bin/skelconv $path$file"_s"$persistence".up.NDskl" -outName $path$file"_s"$persistence".up.NDskl" -smooth $smoothing -to NDskl_ascii


""")
    fichier.write("exit 0")        
    fichier.close()
    os.system(f"sbatch /home/fcastillo/bash/mse_{n}_{i}.sh")
