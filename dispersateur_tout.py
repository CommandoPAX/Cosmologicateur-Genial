import os
import sys
from astropy.io import fits
import numpy as np
import time

pre = "/data100/fcastillo/RESULT/"
snapshots = ["benchM","NG_F500","G_m500","NG_F500_m500","NG_Fminus500","NG_Fminus500_m500", "G_ViVi", "NG_ViVi", "NG_Fminus500_ViVi"]
for n in range(6,9):
    fichier = open(f"../bash/mse_{n}.sh","w")
    fichier.write(f"""#!/bin/bash
#SBATCH --job-name=skl2_{n}
#SBATCH --mem=500gb
#SBATCH --time=48:00:00
#SBATCH --nodelist=i{21+n}
#SBATCH --ntasks=128
#SBATCH --output=/home/fcastillo/logs/skl2_{n}.out
""")
    fichier.write("""
module () {
  eval $(/usr/bin/modulecmd bash $*)
}

source /etc/profile.d/modules.sh
module purge
module load disperse/0.9.24""")

    for i in range(5):

        input_ = pre + snapshots[n]+"/"+str(i)+"_densite"




        fichier.write(f"""
path_dat="{pre+snapshots[n]}/"
path="{pre+snapshots[n]}/"

file="{i}_densite_smooth2"

persistence=0.1
smoothing=1

/softs/disperse/0.9.24/bin/mse $path_dat$file".fits" -nthreads 128 -cut $persistence -upSkl -manifolds -outName $path$file -periodicity 0 -forceLoops
/softs/disperse/0.9.24/bin/skelconv $path$file"_c"$persistence".up.NDskl" -outName $path$file"_c"$persistence".up.NDskl" -smooth $smoothing -to NDskl_ascii

/softs/disperse/0.9.24/bin/mse $path_dat$file".fits" -nthreads 128 -cut $persistence -loadMSC $path$file".MSC" -dumpManifolds J0a -upSkl -manifolds -outName $path$file -periodicity 0 -forceLoops
/softs/disperse/0.9.24/bin/netconv $path$file"_c"$persistence"_manifolds_J0a.NDnet" -outName $path$file"_c"$persistence"_manifolds_J0a.NDnet" -smooth $smoothing -to NDnet_ascii

/softs/disperse/0.9.24/bin/mse $path_dat$file".fits" -nthreads 128 -cut $persistence -loadMSC $path$file".MSC" -dumpManifolds J1a -upSkl -manifolds -outName $path$file -periodicity 0 -forceLoops
/softs/disperse/0.9.24/bin/netconv $path$file"_c"$persistence"_manifolds_J1a.NDnet" -outName $path$file"_c"$persistence"_manifolds_J1a.NDnet" -smooth $smoothing -to NDnet_ascii


""")
    fichier.write("exit 0")        
    fichier.close()
    #os.system(f"sbatch ./bash/mse_{n}_{i}.sh")
