import numpy as np

from pyMIN import*
import matplotlib.pyplot as plt
from math import*
import sys
from astropy.io import fits

n = int(sys.argv[1])
i = int(sys.argv[2])

pre = "/data100/fcastillo/RESULT/"
snapshots = ["benchM","NG_F500","G_m500","NG_F500_m500","NG_Fminus500","NG_Fminus500_m500"]
data = pre + snapshots[n]+"/"+str(i)+"_densite.fits"

hdul = fits.open(data)
data = hdul[0].data
hdul.close()


v0,v1,v2,v3 = calculateMFs(data)


plt.plot(v0,label="v0")
plt.legend()


plt.plot(v1,label="v1")
plt.legend()


plt.plot(v2,label="v2")
plt.legend()


plt.plot(v3,label="v3")

plt.legend()
plt.savefig(pre+snapshots[n]+"/minkowski_"+str(i)+".png")
