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
    
    #snapshots = ["benchM","NG_F500","G_m500","NG_F500_m500","NG_Fminus500","NG_Fminus500_m500"]
    #labels = [r"$\Lambda$CDM", "fnl = -500", "m = 500 eV", "WDM & fnl = -500", "fnl = 500", "WDM & fnl = 500"]

    snapshots = ["benchM","NG_F500","G_ViVi","NG_ViVi","NG_Fminus500","NG_Fminus500_ViVi"]
    labels = [r"$\Lambda$CDM", "fnl = -500", r"$m_{\rm WDM} = 10$ ev, $f_{\rm WDM}$ = 2%", "fnl = -500 & mixed DM", "fnl = 500", "fnl = 500 & mixed DM"]

    lss = ["-", "-", "-.", "--", "-", "--"]
    couleurs = ["blue", "orange", "green", "orange", "fuchsia", "fuchsia"]

    plt.figure(figsize=(14,10))
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
    


    outer = gridspec.GridSpec(nrows=2, ncols=2)

    axs = []
    for row in range(2):
        for col in range(2):
            inner = gridspec.GridSpecFromSubplotSpec(nrows=2, ncols=1, subplot_spec=outer[row, col], hspace=0)
            axs += [plt.subplot(cell) for cell in inner]



    for i in [1,2,3,4]:
        lcdm = np.load(f"/data100/fcastillo/RESULT/{snapshots[0]}/{i}_densite_smooth2_c0.1_connect_fil.txt.npy")
        hist_lcdm = np.histogram(lcdm,  density= True, range = [0, 10], bins=nbins)
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
                        
                    connect = np.load(f"/data100/fcastillo/RESULT/{snapshots[j]}/{i}_densite_smooth2_c0.1_connect_fil.txt.npy")
                    #print(np.shape(data[p]))
                    #print(data)
                    hist = np.histogram(connect, density= True, range = [0, 10], bins=nbins)
                    hist = hist[0]
                    

                    if d == 0 : 
                        axes.plot(hist, color= couleur, ls = ls, label=label)
                        axes.set_ylabel("Probability")
                        axes.xaxis.set_visible(False)
                        axes.set_ylim(0,0.25)
                    else : 
                        axes.plot(hist-hist_lcdm, color=couleurs[j], ls=ls,label=labels[j])
                        axes.set_ylabel(r"$\Delta$")
                    if d ==1 : axes.set_xlabel("Connectivity")

                    if j == 5 and i == 1 and d == 0: 
                        axes.legend() 
                
                except : pass

    if i == 0:
        plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)

    plt.savefig(f"connect_{nbins}.pdf")
    plt.savefig(f"connect_{nbins}.png")

    plt.clf()

    plt.figure()

    for i in range(6):
        moyennes = []
        err = []

        ls = lss[i]
        couleur = couleurs[i]
        label = labels[i]
    
        for j in range(5):
            try :
                    
                connect = np.load(f"/data100/fcastillo/RESULT/{snapshots[i]}/{j}_densite_smooth2_c0.1_connect_fil.txt.npy")
                moyennes.append(np.mean(connect))
                err.append(1/sqrt(len(connect))) *np.std(connect)
                
            except : pass

        print(moyennes)
        axes = plt.gca()

        if len(moyennes )== 5:
            a = 1/(1+np.array([32,3,1,0.25,0]))
        else :
            a = 1/(1+np.array([3,1,0.25,0]))


        #plt.scatter(np.log(np.array([32,3,1,0.25,0][:len(moyennes)])), moyennes, color=couleur)
        plt.errorbar(a, moyennes,ls=ls, color=couleur, label=label, yerr = err)

        axes.set_xlabel(r"$a$")
        axes.set_ylabel("Mean connectivity")
        axes.invert_xaxis()

        plt.legend()
    
    plt.savefig(f"connect_moyenne.pdf")
    plt.savefig(f"connect_moyenne.png")

