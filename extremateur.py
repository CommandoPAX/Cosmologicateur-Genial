from py_extrema.extrema import ExtremaFinder, CriticalPoints
from FyeldGenerator import generate_field
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import ipyvolume.pylab as p34
from astropy.io import fits
from scipy.spatial import KDTree


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

R = 5

try : 
    result = np.load(f"/data100/fcastillo/RESULT/extrema/extrema_{n}_{i}_{R}.txt.npy")
    print("trouve")

except: 
    print("pas trouve")
    hdul = fits.open(data)
    data = hdul[0].data
    hdul.close()
    print(np.shape(data))
    data = data.astype(np.float64)


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
    np.save(f"/data100/fcastillo/RESULT/extrema/extrema_{n}_{i}_{R}.txt", result)
    print("Extrema calcules")


#for j in range(4):
#    for k in range(j, 4):

if True :
    if True :
        j = 3
        k = 3
        print(j, k)
        result_k = result[result[:,3]==k][:,0:3]
        result_j = result[result[:,3]==j][:,0:3]

        data_random = np.load(f"/data100/fcastillo/RESULT/extrema/random.txt.npy")[:len(result_k),:]

        print(result_k)

        tree_k = KDTree(result_k)
        tree_j = KDTree(result_j)
        tree_r = KDTree(data_random)

        print("points charges")

        r_bins = np.linspace(0, 180, 20)  # 20 intervalles entre 0 et 0.2
        print(r_bins)
        r_mid = (r_bins[:-1] + r_bins[1:]) / 2  # Centres des intervalles

        # Compter les paires DD (données-données)

        DD_counts = np.array([len(tree_k.query_ball_tree(tree_j, r)) for r in r_bins[1:]])


        # Compter les paires DR (données-aléatoires)
        DRk_counts = np.array([len(tree_k.query_ball_tree(tree_r, r)) for r in r_bins[1:]])
        DRj_counts = np.array([len(tree_j.query_ball_tree(tree_r, r)) for r in r_bins[1:]])


        # Éviter la division par zéro
        DRk_counts[DRk_counts == 0] = 1
        DRj_counts[DRj_counts == 0] = 1

        # Calcul de ξ(r) selon Davis & Peebles
        correlation = DD_counts / np.sqrt(DRk_counts * DRj_counts) - 1

        print(DD_counts)
        print(DRk_counts)
        print(DRj_counts)
        print(correlation, np.array(correlation))

        np.save(f"/data100/fcastillo/RESULT/extrema/snapshot_{n}_{i}_zeta_{k}_{j}_s{R}.txt", np.array(correlation))
