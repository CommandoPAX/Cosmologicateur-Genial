import numpy as np
import os
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from astropy.io import fits
from math import*
from matplotlib import gridspec

if __name__ == "__main__" :

    Redshifts = {
        0 : 32,
        1 : 3,
        2 : 1,
        3 : 0.25,
        4 : 0
    }
    
    snapshots = ["benchM","NG_F500","G_m500","NG_F500_m500","NG_Fminus500","NG_Fminus500_m500","G_ViVi"]
    labels = [r"$\Lambda$CDM", "fnl = -500", "m = 500 eV", "WDM & fnl = -500", "fnl = 500", "WDM & fnl = 500",r"$m_{\rm WDM} = 10$ ev, $f_{\rm WDM}$ = 2%"]

    #snapshots = ["benchM","NG_F500","G_ViVi","NG_ViVi","NG_Fminus500","NG_Fminus500_ViVi"]
    #labels = [r"$\Lambda$CDM", "fnl = -500", r"$m_{\rm WDM} = 10$ ev, $f_{\rm WDM}$ = 2%", "fnl = -500 & mixed DM", "fnl = 500", "fnl = 500 & mixed DM"]


    ls = ["-", "-", "-.", "--", "-", "--","-"]
    couleurs = ["blue", "orange", "green", "orange", "fuchsia", "fuchsia","green"]

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

    outer = gridspec.GridSpec(nrows=2, ncols=4)

    axs = []
    for row in range(2):
        for col in range(4):
            inner = gridspec.GridSpecFromSubplotSpec(nrows=2, ncols=1, subplot_spec=outer[row, col], hspace=0)
            axs += [plt.subplot(cell) for cell in inner]

    for i in [2,4] :



        for p in range(4):
            lcdm = np.load(f"/home/fcastillo/minkowski/minkowski_{0}_{i}.txt.npy")
            if i in [2,4]:lcdmzoom = np.load(f"/data100/fcastillo/RESULT/benchM/{i}_minkowski_zoom.txt.npy")

            X = np.linspace(-3,3,61)
            if i in [2,4]: 
                X1 = X[X<-1]
                X2 = np.linspace(-1,1,101)
                X3 = X[X>1]

                print(np.shape(X1))

                X = np.concatenate((X1,X2,X3))

            if p == 0 and i in [2,4] : X = X2

            for d in range(2):
                place = places[str(p) + str(d)]

                print(axs)

                axes = axs[(place-1)+(i-2)*4]

                if d == 0 : 
                    axes.title.set_text (r"$z = $"+str(Redshifts[i]))


                for j in range(7):
                    if j in [0,2,6]:
                        lcdm = np.load(f"/home/fcastillo/minkowski/minkowski_{0}_{i}.txt.npy")

                        try:
                            data = np.load(f"/home/fcastillo/minkowski/minkowski_{j}_{i}.txt.npy")
                        except :
                            data= np.load(f"/data100/fcastillo/RESULT/{snapshots[j]}/{i}_minkowski.txt.npy")
                        if i in [2,4]:
                            datazoom = np.load(f"/data100/fcastillo/RESULT/{snapshots[j]}/{i}_minkowski_zoom.txt.npy")

                            data1 = data[:,:20]
                            data2 = data[:,41:]
                            print(np.shape(data1), np.shape(data1),np.shape(data2),np.shape(datazoom))

                            data = np.concatenate([data1,datazoom, data2],axis=1)


                            lcdm = np.concatenate([lcdm[:,:20],lcdmzoom, lcdm[:,41:]],axis=1)

                            if p == 0 : 
                                data = datazoom
                                lcdm = lcdmzoom


                        #print(np.shape(data[p]))
                        #print(data)
                        if d == 0 : 
                            axes.plot(X, data[p], color=couleurs[j], ls=ls[j],label=labels[j])
                            axes.set_ylabel(rf"v$_{p}$")
                            axes.xaxis.set_visible(False)
                        else : 
                            axes.plot(X, data[p] - lcdm[p], color=couleurs[j], ls=ls[j],label=labels[j])
                            axes.set_ylabel(r"$\Delta$")
                        if d ==1 : axes.set_xlabel(r"threshold [$\sigma$]")
                        if j == 6 and i == 0 and d == 0: 
                            axes.legend() 

                        if p == 0 and i in [2,4] : axes.set_xlim(-1,1)
                        else : axes.set_xlim (-3,3)

        if i == 0:
            plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)

        plt.savefig(f"v_dm.pdf")
        plt.savefig(f"v_dm.png")

        #plt.clf()

