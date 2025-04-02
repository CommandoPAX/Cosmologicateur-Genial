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


R = 2

for n in range(6,9):
    for i in range(5):


        print(n, i)

        pre = "/data100/fcastillo/RESULT/"
        snapshots = ["benchM","NG_F500","G_m500","NG_F500_m500","NG_Fminus500","NG_Fminus500_m500", "G_ViVi", "NG_ViVi", "NG_Fminus500_ViVi"]
        data = pre + snapshots[n]+"/"+str(i)+"_densite.fits"
        data0 = pre + snapshots[n]+"/"+str(0)+"_densite.fits"

        #R = 5

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

        ef = ExtremaFinder(field.astype(np.float64))

        field = ef.smooth(R)

        if i != 4 : Npoints = 100
        else : Npoints = 200

        threshold = np.linspace(-6,6,Npoints)

        x = np.arange(512)
        y = np.arange(512)
        z = np.arange(512)

        interpolateur = RegularGridInterpolator((x, y, z), field, bounds_error=False, fill_value=None)

        count = []
        N = 0
        delta = threshold[1] - threshold[0]

        densites = []

        field = field

        for j in range(4):
            type_t = result[result[:, 3] == j]  
            positions = type_t[:, :3]  
            densites_interpolees = interpolateur(positions)
            densites_interpolees = densites_interpolees/np.std(densites_interpolees)


            seuil_haut = np.percentile(densites_interpolees, 95)  # 5% des points les plus hauts
            seuil_bas = np.percentile(densites_interpolees, 5)    # 5% des points les plus bas



            #if j in (0, 1):  # Peaks et filaments : ν > seuil_haut
            #    densites_interpolees = densites_interpolees[densites_interpolees >= seuil_haut]
            #elif j in (2, 3):  # Vides et murs : ν < seuil_bas"""
            #    densites_interpolees = densites_interpolees[(densites_interpolees <= seuil_bas)]
                
            #print(len(type_t), len(densites_interpolees))

            densites.append(densites_interpolees)

            np.save(f"densites_crit_{n}_{i}_s{R}_{j}.txt", np.array(densites_interpolees))


        for t in threshold:
            count_t = []

            for j in range(4):

                densites_interpolees = densites[j]

                points_filtres_t = densites_interpolees[(densites_interpolees> (t-delta)) & (densites_interpolees < t)]

                count_t.append(len(points_filtres_t))
                N += len(points_filtres_t)

            count.append(count_t)
        np.save(f"/data100/fcastillo/RESULT/extrema/snapshot_{n}_{i}_threshold_s{R}.txt", np.array(count)/N)
