from py_extrema.extrema import ExtremaFinder, CriticalPoints
from FyeldGenerator import generate_field
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import ipyvolume.pylab as p34
from astropy.io import fits


from scipy.interpolate import interp1d, RegularGridInterpolator
from collections import defaultdict
from tqdm import tqdm
import pandas as pd


from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt


import sys

n = int(sys.argv[1])
i = int(sys.argv[2])

print(n, i)

pre = "/data100/fcastillo/RESULT/"
snapshots = ["benchM","NG_F500","G_m500","NG_F500_m500","NG_Fminus500","NG_Fminus500_m500"]
data = pre + snapshots[n]+"/"+str(i)+"_densite.fits"

hdul = fits.open(data)
data = hdul[0].data
hdul.close()
print(np.shape(data))
data = data.astype(np.float)

R = 10

ef = ExtremaFinder(data, nthreads=32, loglevel=30)
ef.find_extrema(R)
pos = ef.extrema[R].pos
kind = ef.extrema[R].kind

result = []
for p in range(len(pos)):
    x = float(pos[p][0])
    y = float(pos[p][1])
    z = float(pos[p][2])
    type = kind[p]
    result.append(np.array([x,y,z,type]))


print(np.shape(result))
np.save(f"extrema_{n}_{i}.txt", result)
print("ok")
