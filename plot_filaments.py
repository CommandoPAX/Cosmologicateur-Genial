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

    nbins = 15

    for i in range(1,5):
        plt.subplot(2,2,i)

        axes = plt.gca()

        axes.title.set_text (f"z = {Redshifts[i]}")

        for j in range(6):
                
                ls = lss[j]
                couleur = couleurs[j]
                label = labels[j]

                try :
                    longueurs = np.load(f"/data100/fcastillo/RESULT/{snapshots[j]}/{i}_densite_smooth2_c0.1_len_fil.txt.npy")

                    hist = np.histogram(longueurs, density= True, range = [0, 30], bins=nbins)
                    hist = hist[0]
                    
                    axes.plot(hist, color= couleur, ls = ls, label=label)
                    axes.set_xlabel("length [Mpc / h]")
                    #plt.xscale("log")
                    axes.set_ylabel("PDF")
                    axes.set_ylim(0,0.08)
                    #axes.axvline(np.median(longueurs), color= couleur, ls = ls)
                except:
                     print(j, i)
                    
                if j == 5 and i == 1: plt.legend() 

    plt.savefig(f"len_{nbins}.pdf")

    plt.clf()

    plt.figure()

    for i in range(6):
        moyennes = []
        err = []

        ls = lss[i]
        couleur = couleurs[i]
        label = labels[i]
    
        for j in range(1,5):
            try :
                    
                long = np.load(f"/data100/fcastillo/RESULT/{snapshots[i]}/{j}_densite_smooth2_c0.1_len_fil.txt.npy")
                moyennes.append(np.median(long))
                err.append(1/sqrt(len(long))) *np.std(long)
                
            except : pass

        print(moyennes)
        axes = plt.gca()

        a = 1/(1+np.array([3,1,0.25,0]))


        #plt.scatter(np.log(np.array([32,3,1,0.25,0][:len(moyennes)])), moyennes, color=couleur)
        plt.errorbar(a, moyennes,ls=ls, color=couleur, label=label, yerr = err)

        axes.set_xlabel(r"$a$")
        axes.set_ylabel("Median length")
        axes.invert_xaxis()

        plt.legend()
    
    plt.savefig(f"len_moyenne.pdf")
    plt.savefig(f"len_moyenne.png")

