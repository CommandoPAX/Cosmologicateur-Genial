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
    
    #snapshots = ["benchM","NG_F500","G_m500","NG_F500_m500","NG_Fminus500","NG_Fminus500_m500"]
    #labels = ["LCDM", "fnl = -500", "m = 500 eV", "WDM & fnl = -500", "fnl = 500", "WDM & fnl = 500"]

    snapshots = ["benchM","NG_F500","G_ViVi","NG_ViVi","NG_Fminus500","NG_Fminus500_ViVi"]
    labels = [r"$\Lambda$CDM", "fnl = -500", r"$m_{\rm WDM} = 10$ ev, $f_{\rm WDM}$ = 2%", "fnl = -500 & mixed DM", "fnl = 500", "fnl = 500 & mixed DM"]


    #snapshots = ["benchM", "NEDE","NsPNG_EDE_F1833", "G_ViVi"]
    #labels = ["LCDM", "EDE",  "fnl = -1100 & EDE","mixed DM"]

    lss = ["-", "-", "-.",  "--","-","--"]
    couleurs = ["blue", "orange", "green","orange","fuchsia","fuchsia"]

    #lss = ["-", "-.",  "--","-."]
    #couleurs = ["blue", "red", "#FF9900","green"]
    #ED1C24

    plt.figure(figsize=(14,10))

    nbins = 15

    z= [15,12, 10, 8, 5,3,1,0.5,0.25,0]
    indices_z = [5,6,8,9]

    k = 0

    for i in indices_z:
        k += 1
        plt.subplot(2,2,k)

        axes = plt.gca()

        axes.title.set_text (f"z = {z[i]}")

        for j in range(len(snapshots)):
                
                     
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
                    longueurs = np.load(f"/data100/fcastillo/RESULT/{snapshots[j]}/{k}_densite_smooth2_c0.1_len_fil.txt.npy")

                    hist = np.histogram(longueurs, density= True, range = [0, 30], bins=nbins)
                    hist = hist[0]

                    axes.plot(hist, color= couleur, ls = ls, label=label)
                    axes.set_xlabel("length [Mpc / h]")
                    #plt.xscale("log")
                    axes.set_ylabel("PDF")
                    axes.set_ylim(0,0.08)
                    #axes.axvline(np.median(longueurs), color= couleur, ls = ls)
                    
                if j == 5 and i == 1: plt.legend() 

    plt.savefig(f"len_{nbins}_EDE.pdf")

    plt.clf()

    plt.figure()

    for i in range(len(snapshots)):
        k = 1
        moyennes = []
        err = []

        ls = lss[i]
        couleur = couleurs[i]
        label = labels[i]
    
        for j in indices_z:

            try : long = np.load(f"/data100/fcastillo/RESULT/{snapshots[i]}/{k}_densite_smooth2_c0.1_len_fil.txt.npy")
            except :     long = np.load(f"/data100/fcastillo/RESULT/{snapshots[i]}/{j}_densite_smooth2_c0.1_len_fil.txt.npy")
            moyennes.append(np.mean(long))
            err.append(1/sqrt(len(long)) *np.std(long))

            k += 1
        

        print(moyennes)
        axes = plt.gca()

        a = 1/(1+np.array([3,1,0.25,0]))


        #plt.scatter(np.log(np.array([32,3,1,0.25,0][:len(moyennes)])), moyennes, color=couleur)
        plt.errorbar(np.log10(1+np.array([3,1,0.25,0])), moyennes,ls=ls, color=couleur, label=label, yerr = err)

        axes.set_xlabel(r"$\log (1+z)$")
        axes.set_ylabel("Mean length [Mpc / h]")
        #axes.invert_xaxis()

        plt.legend()
    axes.invert_xaxis()
    plt.savefig(f"len_moyenne_EDE.pdf")
    plt.savefig(f"len_moyenne_EDE.png")

