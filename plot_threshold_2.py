import numpy as np
import os
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from astropy.io import fits
from math import*
from matplotlib import gridspec

snapshots = ["benchM","NG_F500","G_m500","NG_F500_m500","NG_Fminus500","NG_Fminus500_m500","G_ViVi"]
labels = [r"$\Lambda$CDM", "fnl = -500", "m = 500 eV", "WDM & fnl = -500", "fnl = 500", "WDM & fnl = 500",r"$m_{\rm WDM} = 10$ ev, $f_{\rm WDM}$ = 2%"]


if __name__ == "__main__" :

    Redshifts = {
        0 : 32,
        1 : 3,
        2 : 1,
        3 : 0.25,
        4 : 0
    }
    
    #snapshots = ["benchM","NG_F500","G_m500","NG_F500_m500","NG_Fminus500","NG_Fminus500_m500","G_ViVi"]
    #labels = [r"$\Lambda$CDM", "fnl = -500", "m = 500 eV", "WDM & fnl = -500", "fnl = 500", "WDM & fnl = 500",r"$m_{\rm WDM} = 10$ ev, $f_{\rm WDM}$ = 2%"]

    snapshots = ["benchM","NG_F500","G_ViVi","NG_ViVi","NG_Fminus500","NG_Fminus500_ViVi"]
    labels = [r"$\Lambda$CDM", "fnl = -500", r"$m_{\rm WDM} = 10$ ev, $f_{\rm WDM}$ = 2%", "fnl = -500 & mixed DM", "fnl = 500", "fnl = 500 & mixed DM"]


    lss = ["-", "-", "-.", "--", "-", "--","-"]
    couleurs = ["blue", "orange", "green", "orange", "fuchsia", "fuchsia"]

    plt.figure(figsize=(14,10))
    places = {
        "00" : 1,
        "01" : 2,
        "10" : 3,
        "11" : 4,
        "20" : 5,
        "21" : 6,
        "30" : 7,
        "31" : 8
    }
    

    #plt.title(  rf"v$_{p}$")
    plt.axis("off")
    #plt.tight_layout()

    outer = gridspec.GridSpec(nrows=4, ncols=4)
    R = 2

    axs = []
    for row in range(4):
        for col in range(4):
            inner = gridspec.GridSpecFromSubplotSpec(nrows=2, ncols=1, subplot_spec=outer[row, col], hspace=0)
            axs += [plt.subplot(cell) for cell in inner]

    for i in [0,1,2,4] :



        for p in range(4):
            lcdm = np.load(f"/data100/fcastillo/RESULT/extrema/snapshot_{0}_{i}_threshold_s{R}.txt.npy")



            for d in range(2):
                place = places[str(p) + str(d)]

                print(axs)

                axes = axs[(place-1)+(min(i,3))*8]

                if d == 0 : 
                    axes.title.set_text (rf"$v_{p}$,  $z = $"+str(Redshifts[i]))


                for j in range(6):
                    count = np.load(f"/data100/fcastillo/RESULT/extrema/snapshot_{j}_{i}_threshold_s{R}.txt.npy")
                    #zeta[0] = 0
                    Npoints = len(count)
                    X = np.linspace(-6,6,Npoints)

                    couleur = couleurs[j]
                    label = labels[j]
                    ls = lss[j]



                    #axes.plot(X,zeta + 1, color= couleur, ls = ls, label=label)
                    #plt.xscale("log")
                    #axes.set_ylim(0,0.35)

                    if i == 0  and p in (0,1) : axes.set_xlim(-6,6)
                    if i == 1  and p in (0,1) : axes.set_xlim(-6,6)
                    if i == 0  and p in (2,3) : axes.set_xlim(-6,6)
                    if i == 1  and p in (2,3) : axes.set_xlim(-6,6)
                    if i == 2 and p in (0,1) : axes.set_xlim(-2,2)
                    if i == 4 and p in (0,1) : axes.set_xlim(-2,2)
                    if i == 2 and p in (2,3) : axes.set_xlim(-4,0)
                    if i == 4 and p in (2,3) : axes.set_xlim(-4,0)
                    if d == 0 : 
                        axes.plot(X, count[:,p] , color=couleur, ls=ls,label=label)
                        axes.set_ylabel(r"$\frac{1}{N} \frac{dN}{d\nu}$")
                        axes.xaxis.set_visible(False)
                    else : 
                        axes.plot(X, count[:,p]-lcdm[:,p], color=couleur, ls=ls,label=label)
                        axes.set_ylabel(r"$\Delta$")
                    if d ==1 : 
                        axes.set_xlabel(r"$\nu [\sigma]$")
                        #axes.set_xlim(-6,6)

                    if j == 5 and i == 0 and d == 0: 
                        axes.legend() 
                #except: pass

        if i == 0:
            plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)



            
    plt.savefig(f"crit_threshold_s{R}.pdf")
    plt.savefig(f"crit_threshold_s{R}.png")
    
