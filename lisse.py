from astropy.io import fits
import numpy as np
from py_extrema.extrema import ExtremaFinder, CriticalPoints

pre = "/data100/fcastillo/RESULT/"
snapshots = ["SIM0N_BIG","SIM1N_BIG","SIM2N_BIG","SIM3N_BIG"]

for R in [1,2,5]:
    for n in range(4):
        for i in range(4):
            print(R, n, i)
            data = pre + snapshots[n]+"/"+str(i)+"_densite.fits"
            hdul = fits.open(data)
            data = hdul[0].data
            hdul.close()

            ef = ExtremaFinder(data.astype(np.float32), loglevel=30)  
            field = ef.smooth(R*2)

            header = fits.Header()
            header['COMMENT'] = 'Champ de densite'
            header['NAXIS'] = 3  
            header['NAXIS1'] = 1024
            header['NAXIS2'] = 1024 
            header['NAXIS3'] = 1024

            hdu = fits.PrimaryHDU(data=field, header=header)
            hdu.writeto(f"/data100/fcastillo/RESULT/{snapshots[n]}/{i}_densite_smooth{R}.fits", overwrite=True)
