import numpy as np
import os
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from astropy.io import fits
from math import*
from skeletonnateur_hpc import*
from matplotlib import gridspec
import matplotlib


matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams["figure.facecolor"]='w'
matplotlib.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'


if __name__ == "__main__" :

    Redshifts = {
        0 : 32,
        1 : 3,
        2 : 1,
        3 : 0.25,
        4 : 0
    }
    
    #snapshots = ["benchM","NG_F500","G_m500","NG_F500_m500","NG_Fminus500","NG_Fminus500_m500"]
    #labels = [r"$\Lambda$CDM", "fnl = -500", "m = 500 eV", "WDM & fnl = -500", "fnl = 500", "WDM & fnl = 500"]

    #snapshots = ["benchM","NG_F500","G_ViVi","NG_ViVi","NG_Fminus500","NG_Fminus500_ViVi"]
    #labels = [r"$\Lambda$CDM", "fnl = -500", r"$m_{\rm WDM} = 10$ ev, $f_{\rm WDM}$ = 2%", "fnl = -500 & mixed DM", "fnl = 500", "fnl = 500 & mixed DM"]

    #lss = ["-", "-", "-.", "--", "-", "--"]
    #couleurs = ["blue", "orange", "green", "orange", "fuchsia", "fuchsia"]


    #snapshots = ["benchM", "NEDE","NsPNG_EDE_F1833", "G_ViVi"]
    #labels = ["LCDM", "EDE",  "fnl = -1100 & EDE","mixed DM"]


    #lss = ["-", "-.",  "--","-."]
    #couleurs = ["blue", "red", "#FF9900","green"]

    snapshots = ["benchM", "NsPNG_EDE_F500","NsPNG_EDE_F1833", "NG_ViVi","NG_Fminus500_ViVi"]
    labels = [r"$\Lambda$CDM", r"$f_{\rm NL} = -300~\&~{\rm EDE}$",  r"$f_{\rm NL} = -1100~\&~{\rm EDE}$", r"$f_{\rm NL} = -500~\&~{\rm mixed~DM}$", r"$f_{\rm NL} = 500~\&~{\rm mixed~DM}$"]


    lss = ["-", ":", "--",  "--","--"]
    couleurs = ["blue", "darkred", "darkred","darkorange","violet"]

    snapshots = ["benchM", "G_ViVi","NG_F500", "NG_Fminus500","NEDE"]
    labels = [r"$\Lambda$CDM", r"$m_{\rm WDM} = 10  {\rm eV}, f_{\rm WDM} = 2 \%$",  r"$f_{\rm NL} = -500$", r"$f_{\rm NL} = 500$ ", r"${\rm EDE}$"]


    lss = ["-", "-.", "-",  "-","-"]
    couleurs = ["blue", "green", "darkorange","violet","darkred"]



    plt.figure(figsize=(6,6))
    X = np.linspace(-3,3,61)
    places = {
        "10" : 1,
        "11" : 2,
        "20" : 3,
        "21" : 4,
        "30" : 5,
        "31" : 6,
        "40" : 7,
        "41" : 8
    }

    nbins = 10

    z= [15,12, 10, 8, 5,3,1,0.5,0.25,0]
    indices_z = [5,6,8,9]
    

    k = 0
    for i in [4]:
        k += 1
        z_k = 9
        lcdm = np.load(f"/data100/fcastillo/RESULT/{snapshots[0]}/{i}_densite_smooth2_c0.1_percolation.txt.npy")
        

        axes = plt.gca()

        axes.title.set_text (f"z = {Redshifts[i]}")


        for j in range(len(snapshots)):

            ls = lss[j]
            couleur = couleurs[j]
            label = labels[j]


            try :
                    
                try : percole = np.load(f"/data100/fcastillo/RESULT/{snapshots[j]}/{i}_densite_smooth2_c0.1_percolation.txt.npy")
                except : percole = np.load(f"/data100/fcastillo/RESULT/{snapshots[j]}/{z_k}_densite_smooth2_c0.1_percolation.txt.npy")
                #print(np.shape(data[p]))
                #print(data)
                

                axes.plot(np.arange(0.2,1.2,0.01),percole, color= couleur, ls = ls, label=label)
                axes.set_ylabel(r"$S$")
                axes.set_xlabel(r"$bb$")

                axes.legend() 
            
            except : pass

    if i == 0:
        plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)

    plt.savefig(f"percole.pdf")
    plt.savefig(f"percole.png")



