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
    Ne garde que les points critiques dont la densité est entre seuil_min et seuil_max.
    
    :param points_critiques: np.array de forme (N, 4) contenant [X, Y, Z, type]
    :param champ_densite: np.array de forme (512, 512, 512) contenant le champ de densité
    :param seuil_max: float, seuil supérieur de densité
    :param seuil_min: float, seuil inférieur de densité
    :return: np.array filtré des points critiques
    """
    # Extraire les coordonnées des points critiques
    X, Y, Z = points_critiques[:, 0].astype(int), points_critiques[:, 1].astype(int), points_critiques[:, 2].astype(int)
    
    # Récupérer la densité aux positions des points critiques
    densites = champ_densite[X, Y, Z]
    
    # Vérifier la taille de densites et de points_critiques
    assert densites.shape[0] == points_critiques.shape[0], "Mismatch in number of points"
    
    # Filtrer les points dont la densité est dans l'intervalle donné
    indices_valides = (densites <= seuil_max) & (densites > seuil_min)
    
    # Correction : s'assurer que l'indexation fonctionne
    if len(indices_valides.shape) > 1:
        indices_valides = indices_valides.flatten()  # Assurer une forme 1D

    points_filtres = points_critiques[indices_valides]  # Appliquer le masque booléen
    
    return points_filtres



n = int(sys.argv[1])
i = int(sys.argv[2])


print(n, i)

pre = "/data100/fcastillo/RESULT/"
snapshots = ["benchM","NG_F500","G_m500","NG_F500_m500","NG_Fminus500","NG_Fminus500_m500"]
data = pre + snapshots[n]+"/"+str(i)+"_densite.fits"
data0 = pre + snapshots[n]+"/"+str(0)+"_densite.fits"

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

if i == 0 : Npoints = 100//R**2
if i == 1 : Npoints = 100//R**2
if i == 2 : Npoints = 100//R**2
if i == 3 : Npoints = 100//R**2
if i == 4 : Npoints = 100//R**2

threshold = np.linspace(-4,6,Npoints)

x = np.arange(512)
y = np.arange(512)
z = np.arange(512)

interpolateur = RegularGridInterpolator((x, y, z), field, bounds_error=False, fill_value=None)

count = []
N = 0
delta = threshold[1] - threshold[0]

densites = []

for j in range(4):
    type_t = result[result[:, 3] == j]  # Sélection des points critiques du type j
    positions = type_t[:, :3]  # Exclure la colonne "type"
    densites_interpolees = interpolateur(positions)


    # Détermination des seuils spécifiques à chaque type
    seuil_haut = np.percentile(densites_interpolees, 95)  # 5% des points les plus hauts
    seuil_bas = np.percentile(densites_interpolees, 5)    # 5% des points les plus bas


    # Sélection des points en fonction de leur rareté
    if j in (0, 1):  # Peaks et filaments : ν > seuil_haut
        densites_interpolees = densites_interpolees[densites_interpolees >= seuil_haut]
    elif j in (2, 3):  # Vides et murs : ν < seuil_bas"""
        densites_interpolees = densites_interpolees[(densites_interpolees <= seuil_bas)]# | (densites_interpolees >= seuil_haut)]

    print(len(type_t), len(densites_interpolees))

    densites.append(densites_interpolees)

for t in threshold:
    count_t = []

    for j in range(4):
        densites_interpolees = densites[j]

        # Application du filtrage avec le nouveau σ (calculé après sélection)
        sigma = np.std(field)

        points_filtres_t = densites_interpolees[(densites_interpolees> (t-delta)*sigma) & (densites_interpolees < t*sigma)]

        count_t.append(len(points_filtres_t))
        N += len(points_filtres_t)

    count.append(count_t)
np.save(f"/data100/fcastillo/RESULT/extrema/snapshot_{n}_{i}_threshold_s{R}.txt", np.array(count)/N)
