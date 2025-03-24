from astropy.io import fits
import numpy as np
from py_extrema.extrema import ExtremaFinder, CriticalPoints

pre = "/data100/fcastillo/RESULT/"
snapshots = ["benchM","NG_F500","G_m500","NG_F500_m500","NG_Fminus500","NG_Fminus500_m500"]

for R in [1,2,5]:
    for n in range(6):
        for i in range(5):
            data = pre + snapshots[n]+"/"+str(i)+"_densite.fits"
            hdul = fits.open(data)
            data = hdul[0].data
            hdul.close()

            ef = ExtremaFinder(data, nthreads=32, loglevel=30)  
            field = ef.smooth(R)

            header = fits.Header()
            header['COMMENT'] = 'Champ de densite'
            header['NAXIS'] = 3  
            header['NAXIS1'] = 512
            header['NAXIS2'] = 512  
            header['NAXIS3'] = 512 

            hdu = fits.PrimaryHDU(data=field, header=header)
            hdu.writeto(f"/data100/fcastillo/RESULT/{snapshots[n]}/{i}_densite_smooth{R}.fits", overwrite=True)
