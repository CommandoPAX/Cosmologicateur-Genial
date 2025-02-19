import os

for n in range(6):
    for i in range(5):
        fichier = open("bash_"+str(n)+"_"+str(i)+".bash","w")

        fichier.write(f"""#!/bin/bash
#SBATCH --job-name=yt_{n}_{i}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=10:00:00
#SBATCH --output=/home/fcastillo/yt_{n}_{i}.out
#SBATCH --mail-user=fabien.castillo@etu.unistra.fr
#SBATCH --mail-type=END,FAIL 

module purge
module load intelpython

python /home/fcastillo/Cosmologicateur-Genial/Cosmologicateur.py -n {n} -z {i}

exit 0""")

        fichier.close()
        os.system("sbatch bash_"+str(n)+"_"+str(i)+".bash")