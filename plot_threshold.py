import numpy as np
import os
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from astropy.io import fits
from math import*
from skeletonnateur_hpc import*
from matplotlib import gridspec



if __name__ == "__main__" :

    Redshifts = {
        0 : 32,
        1 : 3,
        2 : 1,
        3 : 0.25,
        4 : 0
    }
    
    plt.figure(figsize=(14,10))
    
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

    snapshots = ["benchM","NG_F500","G_m500","NG_F500_m500","NG_Fminus500","NG_Fminus500_m500"]
    labels = ["LCDM", "fnl = -500", "m = 500 eV", "WDM & fnl = -500", "fnl = 500", "WDM & fnl = 500"]

    lss = ["-", "-", "-.", "--", "-", "--"]
    couleurs = ["blue", "orange", "green", "orange", "fuchsia", "fuchsia"]

    nom_corr = ["P", "F", "W", "V"]

    plt.figure(figsize=(14,10))

    nbins = 40
    R = 1
    x0 = 0
    x1 = 40

    #for a in range(4) :
    #    for b in range(a,4):
    for b in range(4) :
        plt.axis("off")
        #plt.tight_layout()

        outer = gridspec.GridSpec(nrows=2, ncols=2)


        axs = []
        for row in range(2):
            for col in range(2):
                inner = gridspec.GridSpecFromSubplotSpec(nrows=2, ncols=1, subplot_spec=outer[row, col], hspace=0)
                axs += [plt.subplot(cell) for cell in inner]


        for i in [0,1,2,4]:
            

            lcdm = np.load(f"/data100/fcastillo/RESULT/extrema/snapshot_{0}_{i}_threshold_s{R}.txt.npy")
            
            for d in range(2):
                place = places[str(i) + str(d)]

                axes = axs[place-1]

                if d == 0 : 
                    axes.title.set_text (f"{nom_corr[b]} z = {Redshifts[i]}")



                for j in range(6):
                        
                        
                    ls = lss[j]
                    couleur = couleurs[j]
                    label = labels[j]

                    if True :#try:
                        count = np.load(f"/data100/fcastillo/RESULT/extrema/snapshot_{j}_{i}_threshold_s{R}.txt.npy")
                        #zeta[0] = 0

                        if i == 0 : Npoints = 40//R**2
                        if i == 1 : Npoints = 40//R**2
                        if i == 2 : Npoints = 40//R**2
                        if i == 3 : Npoints = 40//R**2
                        if i == 4 : Npoints = 40//R**2

                        X = np.linspace(-4,6,Npoints)


                        #axes.plot(X,zeta + 1, color= couleur, ls = ls, label=label)
                        #plt.xscale("log")
                        #axes.set_ylim(0,0.35)

                        """if i == 0  and b in (0,1) : axes.set_xlim(-4,4)
                        if i == 1  and b in (0,1) : axes.set_xlim(-2,4)
                        if i == 0  and b in (2,3) : axes.set_xlim(-4,2)
                        if i == 1  and b in (2,3) : axes.set_xlim(-2,2)
                        if i == 2 and b in (0,1) : axes.set_xlim(-1,1)
                        if i == 4 and b in (0,1) : axes.set_xlim(-1,1)
                        if i == 2 and b in (2,3) : axes.set_xlim(-0.5,0.5)
                        if i == 4 and b in (2,3) : axes.set_xlim(-0.5,0.5)"""
                        if d == 0 : 
                            axes.plot(X, count[:,b] , color=couleur, ls=ls,label=label)
                            axes.set_ylabel(r"$\frac{1}{N} \frac{dN}{d\nu}$")
                            axes.xaxis.set_visible(False)
                        else : 
                            axes.plot(X, count[:,b]-lcdm[:,b], color=couleur, ls=ls,label=label)
                            axes.set_ylabel(r"$\Delta$")
                        if d ==1 : 
                            axes.set_xlabel(r"$\nu [\sigma]$")
                            axes.set_xlim(-4,6)

                        if j == 5 and i == 0 and d == 0: 
                            axes.legend() 
                    #except: pass

            if i == 0:
                plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)



            
        plt.savefig(f"crit_threshold_{b}_s{R}.pdf")
        plt.savefig(f"crit_threshold_{b}_s{R}.png")
        plt.clf()


    for i in [0,1,2,4]:
        plt.subplot(2,2,min(i+1,4))

        axes = plt.gca()

        axes.title.set_text (f"z = {Redshifts[i]}")

        for j in range(4):
                
                ls = "-"#lss[j]
                couleur = ["orange","green","blue","red"][3-j]
                label = nom_corr[j]#labels[j]

                try :
                    count = np.load(f"/data100/fcastillo/RESULT/extrema/snapshot_{0}_{i}_threshold_s{R}.txt.npy")

                    print(count)

                    axes.plot(X, count[:,j] , color=couleur, ls=ls,label=label)
                    axes.set_xlabel(r"$\nu [\sigma]$")
                    axes.set_ylabel(r"$\frac{1}{N} \frac{dN}{d\nu}$")
                except:
                     print(j, i)
                    
                if j == 3 and i == 0: plt.legend() 

    plt.savefig("extrema_lcdm.pdf")