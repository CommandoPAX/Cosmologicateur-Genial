import os
import sys
from astropy.io import fits
import numpy as np

pre = "../../../data100/fcastillo/RESULT/"
snapshots = ["benchM","NG_F500","G_m500","NG_F500_m500","NG_Fminus500","NG_Fminus500_m500"]

n = int (sys.argv[1])
i = int(sys.argv[2])

fichier = open("mse_"+str(n)+"_"+str(i)+".bash","w")

input_ = pre + snapshots[n]+"/"+str(i)+"_densite"

hdul = fits.open(input_+".fits")
data = hdul[0].data
slice = np.sum(data[:,:,0:2], axis=2)
hdul.close()

header = fits.Header()
header['COMMENT'] = 'Slice'
header['NAXIS'] = 3  
header['NAXIS1'] = data.shape[0]  
header['NAXIS2'] = data.shape[1]  
header['NAXIS3'] = data.shape[2] 

hdu = fits.PrimaryHDU(data=slice, header=header)

hdu.writeto(input_+"_slice.fits", overwrite=True)


fichier.write(f"""#!/bin/bash
#SBATCH --job-name=disperse_{n}_{i}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100gb
#SBATCH --time=10:00:00
#SBATCH --output=/home/fcastillo/disperse_slice_{n}_{i}.out   
#SBATCH --mail-user=fabien.castillo@etu.unistra.fr
#SBATCH --mail-type=ALL 

module purge
module load disperse/0.9.24

/softs/disperse/0.9.24/bin/mse {input_}_slice.fits -cut 3 -upSkl -manifolds -outName {pre+snapshots[n]+"/"+str(i)+"_densite_slice.fits"}
/softs/disperse/0.9.24/bin/skelconv {pre+snapshots[n]+"/"+str(i)+"_densite_slice.fits_c3.up.NDskl"} -smooth 1 -outName {pre+snapshots[n]+"/"+str(i)+"_densite_slice.fits_c3.up.NDskl"} -to NDskl_ascii

exit 0""")

fichier.close()
os.system("sbatch mse_"+str(n)+"_"+str(i)+".bash")