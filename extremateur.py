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

n = sys.argv[1]
i = sys.argv[2]

print(n, i)

pre = "/data100/fcastillo/RESULT/"
snapshots = ["benchM","NG_F500","G_m500","NG_F500_m500","NG_Fminus500","NG_Fminus500_m500"]
data = pre + snapshots[n]+"/"+str(i)+"_densite.fits"

hdul = fits.open(data)
data = hdul[0].data
hdul.close()
print(np.shape(data))

ef = ExtremaFinder(data, nthreads=32, loglevel=30)
ef.find_extrema(10)
curvature = ef.curvature


v0,v1,v2,v3 = calculateMFs(data)

result = np.array([v0,v1,v2,v3])
print(np.shape(result))
np.save(f"minkowski_{n}_{i}.txt", result)
print("ok")
