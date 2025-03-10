import numpy as np
import os
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from astropy.io import fits
from math import*
from skeletonnateur_hpc import*



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

    lss = ["-", "-", "-.", "--", "-", "--"]
    couleurs = ["blue", "orange", "green", "orange", "fuchsia", "fuchsia"]

    plt.figure(figsize=(14,10))

    nbins = 20
    R = 5
    x0 = 0
    x1 = 180

    #for a in range(4) :
    #    for b in range(a,4):
    if True :
        if True :
            a = 3
            b = 3

            for i in [0,1,2,4]:
                plt.subplot(2,2,min(i+1,4))

                axes = plt.gca()

                axes.title.set_text (f"z = {Redshifts[i]}")

                for j in range(6):
                        
                            
                        ls = lss[j]
                        couleur = couleurs[j]
                        label = labels[j]

                        zeta = np.load(f"/data100/fcastillo/RESULT/extrema/snapshot_{j}_{i}_zeta_{a}_{b}_s{R}.txt.npy")
                        print(zeta)

                        X = np.linspace(x0, x1,nbins-1)

                        axes.plot(X,zeta, color= couleur, ls = ls, label=label)
                        axes.set_xlabel("r [Mpc / h]")
                        #plt.xscale("log")
                        axes.set_ylabel(r"$\zeta (r)$")
                        #axes.set_ylim(0,0.35)

                            
                        if j == 5 and i == 1: plt.legend() 


            plt.savefig(f"corr_{a}_{b}_s{R}_nbins{nbins}.pdf")
            plt.savefig(f"corr_{a}_{b}_s{R}_nbins{nbins}.png")
            plt.clf()

