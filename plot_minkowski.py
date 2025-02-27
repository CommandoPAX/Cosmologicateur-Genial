import numpy as np
import os
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from astropy.io import fits
from math import*

if __name__ == "__main__" :

    Redshifts = {
        0 : 32,
        1 : 3,
        2 : 1,
        3 : 0.25,
        4 : 0
    }
    
    snapshots = ["benchM","NG_F500","G_m500","NG_F500_m500","NG_Fminus500","NG_Fminus500_m500"]
    labels = ["LCDM", "fnl = -500", "m = 500 eV", "WDM & fnl = -500", "fnl = 500", "WDM & fnl = 500"]

    ls = ["-", "-", "-.", "--", "-", "--"]
    couleurs = ["blue", "orange", "green", "orange", "fuchsia", "fuchsia"]

    plt.figure(figsize=(14,10))

    for p in range(4) :
        for i in range(1,5):
            plt.title(f"v_{p}")
            plt.subplot(2,2,i)

            axes = plt.gca()

            axes.title.set_text (f"z = {Redshifts[i]}")

            for j in range(6):

                data = np.load(f"minkowski_{i}_{j}.txt.npy")
                axes.plot(data[p], color=couleurs[j], ls=ls[j],label=labels[j])
                if j == 5 and i == 1: plt.legend() 

        plt.savefig(f"v_{p}.pdf")
        plt.clf()

