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


def filtrer_points_critiques(points_critiques, champ_densite, seuil_max, seuil_min):
    """
    Ne garde que les points critiques dont la densité est supérieure à un seuil donné.
    
    :param points_critiques: np.array de forme (N, 4) contenant [X, Y, Z, type]
    :param champ_densite: np.array de forme (512, 512, 512) contenant le champ de densité
    :param seuil: float, seuil de densité
    :return: np.array filtré des points critiques
    """
    # Extraire les coordonnées des points critiques
    X, Y, Z = points_critiques[:, 0].astype(int), points_critiques[:, 1].astype(int), points_critiques[:, 2].astype(int)
    
    # Récupérer la densité aux positions des points critiques
    densites = champ_densite[X, Y, Z]
    
    # Filtrer les points dont la densité dépasse le seuil
    indices_valides = (densites >= seuil_max) & (densites < seuil_min)
    points_filtres = points_critiques[indices_valides]
    
    return points_filtres






n = int(sys.argv[1])
i = int(sys.argv[2])


print(n, i)

pre = "/data100/fcastillo/RESULT/"
snapshots = ["benchM","NG_F500","G_m500","NG_F500_m500","NG_Fminus500","NG_Fminus500_m500"]
data = pre + snapshots[n]+"/"+str(i)+"_densite.fits"

R = 1

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


hdul = fits.open(data)
field = hdul[0].data
hdul.close()

threshold = np.linspace(-4,6,100)
count = []
sigma = np.std(field)
delta = threshold[1] - threshold[0]

for t in threshold :
    points_filtres = filtrer_points_critiques(result, field, t*sigma, (t-delta)*sigma)

    count_t = []

    for j in range(4) :

        count_t.append(len(points_filtres[points_filtres[:,3]==j]))
    
    count.append(count_t)



np.save(f"/data100/fcastillo/RESULT/extrema/snapshot_{n}_{i}_threshold_s{R}.txt", np.array(count))
