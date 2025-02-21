import numpy as np
from astropy.io import fits
import sys

n = sys.argv[1]
i = sys.argv[2]

pre = "../../../data100/fcastillo/RESULT/"
snapshots = ["benchM","NG_F500","G_m500","NG_F500_m500","NG_Fminus500","NG_Fminus500_m500"]
data = pre + snapshots[n]+"/"+str(i)+"_densite.fits"

hdul = fits.open(data)
delta = hdul[0].data
hdul.close()

n = np.shape(delta)[0]

header = fits.Header()
header['COMMENT'] = 'Champ de densite'
header['NAXIS'] = 3  
header['NAXIS1'] = delta.shape[1]/2  
header['NAXIS2'] = delta.shape[0] /2
header['NAXIS3'] = delta.shape[2] /2

for j in range(8):
 
    if j == 0 :_data = delta[0:n//2, 0:n//2, 0:n//2]
    if j == 1 :_data = delta[0:n//2, 0:n//2, n//2:n]
    if j == 2 :_data = delta[0:n//2, n//2:n, 0:n//2]
    if j == 3 :_data = delta[0:n//2, n//2:n, n//2:n]
    if j == 4 :_data = delta[n//2:n, 0:n//2, 0:n//2]
    if j == 5 :_data = delta[n//2:n, 0:n//2, n//2:n]
    if j == 6 :_data = delta[n//2:n, n//2:n, 0:n//2]
    if j == 7 :_data = delta[n//2:n, n//2:n, n//2:n]


    hdu = fits.PrimaryHDU(data=delta, header=header)

    hdu.writeto(pre + snapshots[n]+"/"+str(i)+"_densite_"+str(j)+".fits", overwrite=True)
   