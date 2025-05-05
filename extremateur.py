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

def count_pairs_query_ball_point(points_A, tree_B, r):
    neighbors_list = tree_B.query_ball_point(points_A, r, workers=-1)
    count = sum(len(n) for n in neighbors_list)

    neighbors_list0 = tree_B.query_ball_point(points_A, 0.1, workers=-1)
    count0 = sum(len(n) for n in neighbors_list0)

    return count - count0

n = int(sys.argv[1])
i = int(sys.argv[2])
P = 5

print(n, i)

pre = "/data100/fcastillo/RESULT/"
snapshots = ["benchM","NG_F500","G_m500","NG_F500_m500","NG_Fminus500","NG_Fminus500_m500", "G_ViVi","NG_ViVi" , "NG_Fminus500_ViVi","NEDE", "NsPNG_EDE_F500","NsPNG_EDE_F1833"]
data = pre + snapshots[n]+"/"+str(i)+"_densite.fits"

R = 2

try : 
    result = np.load(f"/data100/fcastillo/RESULT/extrema/extrema_{n}_{i}_{R}.txt.npy")
    print("trouve")

except: 
    print("pas trouve")
    hdul = fits.open(data)
    data = hdul[0].data
    hdul.close()
    print(np.shape(data))
    data = data.astype(np.float32)


    ef = ExtremaFinder(data, loglevel=30)
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
    result = np.load(f"/data100/fcastillo/RESULT/extrema/extrema_{n}_{i}_{R}.txt.npy")



#for j in range(4):
#    for k in range(j, 4):
data = pre + snapshots[n]+"/"+str(i)+"_densite_smooth2.fits"

hdul = fits.open(data)
field = hdul[0].data
hdul.close()

x = np.arange(512)
y = np.arange(512)
z = np.arange(512)

for j in range(4) :
    for k in range(j,4) :

        print(j, k)
        result_k = result[result[:,3]==k][:,0:3] % 512
        result_j = result[result[:,3]==j][:,0:3] % 512

        Xk, Yk, Zk = result_k[:, 0].astype(int), result_k[:, 1].astype(int), result_k[:, 2].astype(int)
        Xj, Yj, Zj = result_j[:, 0].astype(int), result_j[:, 1].astype(int), result_j[:, 2].astype(int)


        seuil_haut_k = np.percentile(field[Xk,Yk,Zk], 100-P)
        seuil_bas_k = np.percentile(field[Xk,Yk,Zk], P)

        seuil_haut_j = np.percentile(field[Xj,Yj,Zj], 100-P)
        seuil_bas_j = np.percentile(field[Xj,Yj,Zj], P)

        indices_k = field[Xk, Yk, Zk] >= seuil_haut_k
        indices_j = field[Xj, Yj, Zj] >= seuil_haut_j


        if j in (0, 1):  # Peaks et filaments 
            indices_j = field[Xj, Yj, Zj] >= seuil_haut_j
        elif j in (2, 3):  # Vides et murs 
            indices_j = field[Xj, Yj, Zj] <= seuil_bas_j

        if k in (0, 1):  # Peaks et filaments 
            indices_k = field[Xk, Yk, Zk] >= seuil_haut_k
        elif k in (2, 3):  # Vides et murs 
            indices_k = field[Xk, Yk, Zk] <= seuil_bas_k

        data_random = np.random.uniform(low=0,high=512,size=(len(result_k[indices_k])*100,3))

        points_k = result_k[indices_k]
        points_j = result_j[indices_j]

        tree_r = KDTree(data_random, boxsize=512)

        r_small = np.linspace(0.1, 10, 80)  # 10 points entre 0 et 1
        r_large = np.geomspace(10, 20, 20)  # 30 points entre 1 et 40 (logarithmique)
        r_bins = np.concatenate((r_small, r_large))

        DD_counts = np.array([count_pairs_query_ball_point(points_k, KDTree(points_j, boxsize=512), r)
                            for r in r_bins[1:]])

        DRk_counts = np.array([count_pairs_query_ball_point(points_k, tree_r, r)
                            for r in r_bins[1:]])

        DRj_counts = np.array([count_pairs_query_ball_point(points_j, tree_r, r)
                            for r in r_bins[1:]])
        DRk_counts[DRk_counts == 0] = 1
        DRj_counts[DRj_counts == 0] = 1

        correlation = DD_counts / np.sqrt(DRk_counts * DRj_counts) *sqrt(len(data_random) * len(data_random) / (len(result_k[indices_k]) * len(result_j[indices_j]))) - 1

        print(DD_counts)
        print(DRk_counts)
        print(DRj_counts)
        print(correlation, np.array(correlation))

        np.save(f"/data100/fcastillo/RESULT/extrema/snapshot_{n}_{i}_zeta_{k}_{j}_s{R}_P{P}_tres_grand.txt", np.array(correlation))
