from astropy.io import fits
import numpy as np
from py_extrema.extrema import ExtremaFinder, CriticalPoints

pre = "/data100/fcastillo/RESULT/"
snapshots = ["NEDE","NsPNG_F500","NsPNG_F1000","NsPNG_F1833","NsPNG_EDE_F500","NsPNG_EDE_F1000","NsPNG_EDE_F1833"]

for R in [1,2,5]:
    for n in range(7):
        for i in range(10):
            print(R, n, i)
            data = pre + snapshots[n]+"/"+str(i)+"_densite.fits"
            hdul = fits.open(data)
            data = hdul[0].data
            hdul.close()

            ef = ExtremaFinder(data.astype(np.float32), nthreads=32, loglevel=30)  
            field = ef.smooth(R)

            header = fits.Header()
            header['COMMENT'] = 'Champ de densite'
            header['NAXIS'] = 3  
            header['NAXIS1'] = 512
            header['NAXIS2'] = 512  
            header['NAXIS3'] = 512 

            hdu = fits.PrimaryHDU(data=field, header=header)
            hdu.writeto(f"/data100/fcastillo/RESULT/{snapshots[n]}/{i}_densite_smooth{R}.fits", overwrite=True)
