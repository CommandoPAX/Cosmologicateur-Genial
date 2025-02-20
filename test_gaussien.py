import numpy as np
from astropy.io import fits


delta = np.random.normal(size=(64,64,64))
delta /= np.mean(delta, dtype=np.float64);  delta -= 1.0


header = fits.Header()
header['COMMENT'] = 'Champ de densite'
header['NAXIS'] = 3  
header['NAXIS1'] = delta.shape[1]  
header['NAXIS2'] = delta.shape[0]  
header['NAXIS3'] = delta.shape[2] 

hdu = fits.PrimaryHDU(data=delta, header=header)

hdu.writeto(f"gaussien_densite.fits", overwrite=True)
