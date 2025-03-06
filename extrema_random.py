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


data = np.random.uniform(size=(512,512,512))

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
np.save(f"/data100/fcastillo/RESULT/extrema/extrema_random_{R}.txt", result)
print("ok")
