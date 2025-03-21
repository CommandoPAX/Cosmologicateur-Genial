import numpy as np
import os
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from astropy.io import fits
from math import*
from NDnet import*
from matplotlib import gridspec


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

    lss = ["-", "-", "-.", "--", "-", "--"]
    couleurs = ["blue", "orange", "green", "orange", "fuchsia", "fuchsia"]

    plt.figure(figsize=(14,10))
    X = np.linspace(-3,3,61)
    places = {
        "00" : 1,
        "01" : 2,
        "10" : 3,
        "11" : 4,
        "20" : 5,
        "21" : 6,
        "40" : 7,
        "41" : 8
    }

    nbins = 10
    


    outer = gridspec.GridSpec(nrows=2, ncols=2)

    axs = []
    for row in range(2):
        for col in range(2):
            inner = gridspec.GridSpecFromSubplotSpec(nrows=2, ncols=1, subplot_spec=outer[row, col], hspace=0)
            axs += [plt.subplot(cell) for cell in inner]



    for i in [0,1,2,4]:
        lcdm = np.load(f"/data100/fcastillo/RESULT/{snapshots[0]}/{i}_densite_0_c0.1_surf_murs.txt.npy")
        hist_lcdm = np.histogram(lcdm,  density= True, range = [0, 4], bins=nbins)
        hist_lcdm = hist_lcdm[0]


        for d in range(2):
            place = places[str(i) + str(d)]

            print(axs)

            axes = axs[place-1]

            if d == 0 : 
                axes.title.set_text (f"z = {Redshifts[i]}")


            for j in range(6):

                ls = lss[j]
                couleur = couleurs[j]
                label = labels[j]


                try :
                        
                    connect = np.load(f"/data100/fcastillo/RESULT/{snapshots[j]}/{i}_densite_0_c0.1_surf_murs.txt.npy")
                    #print(np.shape(data[p]))
                    #print(data)
                    hist = np.histogram(connect, density= True, range = [0, 10], bins=nbins)
                    hist = hist[0]
                    

                    if d == 0 : 
                        axes.plot(hist, color= couleur, ls = ls, label=label)
                        axes.set_ylabel("Probabilite")
                        axes.xaxis.set_visible(False)
                        axes.set_ylim(0,0.35)
                        axes.set_xlim(0,6)
                    else : 
                        axes.plot(hist-hist_lcdm, color=couleurs[j], ls=ls,label=labels[j])
                        axes.set_ylabel(r"$\Delta$")
                    if d ==1 : axes.set_xlabel("Surface")

                    if j == 5 and i == 0 and d == 0: 
                        axes.legend() 
                
                except : pass

    if i == 0:
        plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)

    plt.savefig(f"Surf_{nbins}.pdf")
    plt.savefig(f"Surf_{nbins}.png")




