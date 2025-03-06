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
from math import sqrt

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
data = data.astype(np.float64)

R = 5

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
np.save(f"extrema/extrema_{n}_{i}_{R}.txt", result)
print("Extrema calcules")

data_random = np.load(f"data100/fcastillo/RESULT/extrema/extrema_random_{R}.txt.npy")

for j in range(4):
    for k in range(4):
        print(j, k)
        result_k = []
        for p in result :
            if p[3] == k : result_k.append(p)


        result_j = []
        for p in result :
            if p[3] == j : result_j.append(p)

        random_k = []
        for p in data_random :
            if p[3] == k : random_k.append(p)


        random_j = []
        for p in data_random :
            if p[3] == j : random_j.append(p)

        Ckj = []
        Rkj = []
        Rjk = []

        for p in result_k :
            x0 = p[0]
            y0 = p[1]
            z0 = p[2]
            for q in random_j :
                x1 = q[0]
                y1 = q[1]
                z1 = q[2]
                Rkj.append(sqrt((x0-x1)**2+(y0-y1)**2 +(z0-z1)**2))       
            for q in result_j :
                x1 = q[0]
                y1 = q[1]
                z1 = q[2]
                Ckj.append(sqrt((x0-x1)**2+(y0-y1)**2 +(z0-z1)**2))   
        for p in random_k :
            x0 = p[0]
            y0 = p[1]
            z0 = p[2]
            for q in result_j :
                x1 = q[0]
                y1 = q[1]
                z1 = q[2]
                Rjk.append(sqrt((x0-x1)**2+(y0-y1)**2 +(z0-z1)**2))       

        np.save(f"/data100/fcastillo/RESULT/extrema/snapshot_{n}_{i}_C_{k}_{j}_s{R}.txt", np.array(Ckj))
        np.save(f"/data100/fcastillo/RESULT/extrema/snapshot_{n}_{i}_R_{k}_{j}_s{R}.txt", np.array(Rkj))
        np.save(f"/data100/fcastillo/RESULT/extrema/snapshot_{n}_{i}_R_{j}_{k}_s{R}.txt", np.array(Rjk))
