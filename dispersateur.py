import os
import sys


pre = "../../../data100/fcastillo/RESULT/"
snapshots = ["benchM","NG_F500","G_m500","NG_F500_m500","NG_Fminus500","NG_Fminus500_m500"]

n = int (sys.argv[1])
i = int(sys.argv[2])

fichier = open("mse_"+str(n)+"_"+str(i)+".bash","w")

input_ = pre + snapshots[n]+"/"+str(i)+"_densite.fits"

fichier.write(f"""#!/bin/bash
#SBATCH --job-name=disperse_{n}_{i}
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --ntasks-per-node=32
#SBATCH --mem=500gb
#SBATCH --time=10:00:00
#SBATCH --output=/home/fcastillo/disperse_{n}_{i}.out   
#SBATCH --mail-user=fabien.castillo@etu.unistra.fr
#SBATCH --mail-type=ALL 

module purge
module load disperse/0.9.24

/softs/disperse/0.9.24/bin/mse {input_} -cut 3 -upSkl -manifolds -nthreads 32 -outName {pre+snapshots[n]+"/"+str(i)+"_densite"}

exit 0""")

fichier.close()
os.system("sbatch mse_"+str(n)+"_"+str(i)+".bash")