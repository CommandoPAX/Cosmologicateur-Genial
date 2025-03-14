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

    nbins = 10

    for i in range(1,5):
        plt.subplot(2,2,i)

        axes = plt.gca()

        axes.title.set_text (f"z = {Redshifts[i]}")

        for j in range(6):
                
                ls = lss[j]
                couleur = couleurs[j]
                label = labels[j]

                try :
                    longueurs = np.load(f"/data100/fcastillo/RESULT/{snapshots[j]}/{i}_densite_0_c0.1_len_fil.txt.npy")

                    hist = np.histogram(longueurs, density= True, range = [0, 10], bins=nbins)
                    hist = hist[0]
                    
                    axes.plot(hist, color= couleur, ls = ls, label=label)
                    axes.set_xlabel("log longueur [Mpc / h]")
                    plt.xscale("log")
                    axes.set_ylabel("Probabilite")
                    axes.set_ylim(0,0.35)
                    axes.axvline(np.median(longueurs), color= couleur, ls = ls)
                except:
                     print(j, i)
                    
                if j == 5 and i == 1: plt.legend() 

    plt.savefig(f"len_{nbins}.pdf")

