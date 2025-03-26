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

def count_pairs(tree_A, tree_B, r):
    global R

    count = 0
    count0 = 0

    if True : #R >= 2 :
        
        count = np.sum(len(neighbors) for neighbors in tree_A.query_ball_tree(tree_B, r))

        count0 = np.sum(len(neighbors) for neighbors in tree_A.query_ball_tree(tree_B, 0.1))

    else :
        for i in range(0, len(tree_A.data), 10000):
            sub_tree_A = KDTree(tree_A.data[i : min(i + 10000,len(tree_A.data))], boxsize=512)
            sub_count = np.sum(len(neighbors) for neighbors in sub_tree_A.query_ball_tree(tree_B, r))
            sub_count_0 = np.sum(len(neighbors) for neighbors in sub_tree_A.query_ball_tree(tree_B, 0.1))
            
            count += sub_count
            count0 += sub_count_0

    return count - count0

n = int(sys.argv[1])
i = int(sys.argv[2])
P = 20

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
    result = np.load(f"/data100/fcastillo/RESULT/extrema/extrema_{n}_{i}_{R}.txt.npy")



#for j in range(4):
#    for k in range(j, 4):

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


        seuil_haut_k = np.percentile(field[Xk,Yk,Zk], P)
        seuil_bas_k = np.percentile(field[Xk,Yk,Zk], 100-P)

        seuil_haut_j = np.percentile(field[Xj,Yj,Zj], P)
        seuil_bas_j = np.percentile(field[Xj,Yj,Zj], 100-P)

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

        data_random = np.load(f"/data100/fcastillo/RESULT/extrema/random.txt.npy")[:len(result_k[indices_k]),:] % 512

        tree_k = KDTree(result_k[indices_k],boxsize=512) 
        tree_j = KDTree(result_j[indices_j],boxsize=512)
        tree_r = KDTree(data_random,boxsize=512)

        print("points charges")

        r_small = np.linspace(0, 5, 80)
        r_large = np.geomspace(5, 40, 20)
        r_bins = np.concatenate((r_small, r_large))


        DD_counts = np.array([count_pairs(tree_k, tree_j, r) for r in r_bins[1:]])

        DRk_counts = np.array([count_pairs(tree_k, tree_r, r) for r in r_bins[1:]])
        DRj_counts = np.array([count_pairs(tree_j, tree_r, r) for r in r_bins[1:]])

        DRk_counts[DRk_counts == 0] = 1
        DRj_counts[DRj_counts == 0] = 1

        correlation = DD_counts / np.sqrt(DRk_counts * DRj_counts) *sqrt(len(data_random) * len(data_random) / (len(result_k) * len(result_j))) - 1

        print(DD_counts)
        print(DRk_counts)
        print(DRj_counts)
        print(correlation, np.array(correlation))

        np.save(f"/data100/fcastillo/RESULT/extrema/snapshot_{n}_{i}_zeta_{k}_{j}_s{R}_P{P}.txt", np.array(correlation))
