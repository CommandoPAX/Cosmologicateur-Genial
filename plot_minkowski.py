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
    X = np.linspace(-3,3,61)
    places = {
        "00" : 1,
        "01" : 3,
        "10" : 2,
        "11" : 4,
        "20" : 5,
        "21" : 7,
        "40" : 6,
        "41" : 8
    }
    


    for p in range(4) :
        plt.title(  f"v_{p}")
        plt.axis("off")
        #plt.tight_layout()
        for i in [0,1,2,4]:
            lcdm = np.load(f"minkowski_{0}_{i}.txt.npy")

            for d in range(2):
                place = places[str(i) + str(d)]
                plt.subplot(4,2,place)

                axes = plt.gca()

                axes.title.set_text (f"z = {Redshifts[i]}")

                for j in range(6):

                    data = np.load(f"minkowski_{j}_{i}.txt.npy")
                    #print(np.shape(data[p]))
                    #print(data)
                    if d == 0 : 
                        axes.plot(X, data[p], color=couleurs[j], ls=ls[j],label=labels[j])
                        axes.set_ylabel(rf"v$_{p}$")
                    else : 
                        axes.plot(X, data[p] - lcdm[p], color=couleurs[j], ls=ls[j],label=labels[j])
                        axes.set_ylabel(r"\Delta")
                    axes.set_xlabel(r"threshold [$\sigma$]")
                    if j == 5 and i == 0: 
                        axes.legend() 

        if i == 0:
            plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)

        plt.savefig(f"v_{p}.pdf")
        plt.savefig(f"v_{p}.png")

        plt.clf()

