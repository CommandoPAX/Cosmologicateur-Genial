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
    labels = [r"$\Lambda$CDM", "fnl = -500", "m = 500 eV", "WDM & fnl = -500", "fnl = 500", "WDM & fnl = 500"]

    ls = ["-", "-", "-.", "--", "-", "--"]
    couleurs = ["blue", "orange", "green", "orange", "fuchsia", "fuchsia"]

    plt.figure(figsize=(14,10))
    X = np.arange(-3,3,61)

    for p in range(4) :
        plt.title(  f"v_{p}")
        #plt.tight_layout()
        for i in [0,1,2,4]:
            plt.subplot(2,2,min(i+1,4))

            axes = plt.gca()

            axes.title.set_text (f"z = {Redshifts[i]}")

            for j in range(6):


                data = np.load(f"minkowski_{j}_{i}.txt.npy")
                print(np.shape(data[p]))
                print(data)
                axes.plot(X, data[p], color=couleurs[j], ls=ls[j],label=labels[j])
                axes.set_xlabel(r"threshold [$\sigma$]")
                if j == 5 and i == 1: 
                    axes.legend() 

        if i == 0:
            plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)

        plt.savefig(f"v_{p}.pdf")
        plt.clf()

