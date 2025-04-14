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
    #labels = ["LCDM", "fnl = -500", "m = 500 eV", "WDM & fnl = -500", "fnl = 500", "WDM & fnl = 500"]

    snapshots = ["benchM","NG_F500","G_ViVi","NG_ViVi","NG_Fminus500","NG_Fminus500_ViVi"]
    labels = [r"$\Lambda$CDM", "fnl = -500", r"$m_{\rm WDM} = 10$ ev, $f_{\rm WDM}$ = 2%", "fnl = -500 & mixed DM", "fnl = 500", "fnl = 500 & mixed DM"]


    #snapshots = ["benchM", "NEDE","NsPNG_EDE_F1833", "G_ViVi"]
    #labels = ["LCDM", "EDE",  "fnl = -1100 & EDE","mixed DM"]

    #lss = ["-", "-.",  "--","-."]
    #couleurs = ["blue", "red", "#FF9900","green"]
    #ED1C24

    snapshots = ["benchM", "G_ViVi","NG_F500", "NG_Fminus500","NEDE"]
    labels = [r"$\Lambda$CDM", r"$m_{\rm WDM} = 10  {\rm eV}, f_{\rm WDM} = 2 \%$",  r"$f_{\rm NL} = -500$", r"$f_{\rm NL} = 500$ ", r"${\rm EDE}$"]


    lss = ["-", "-.", "-",  "-","-"]
    couleurs = ["blue", "green", "darkorange","violet","darkred"]


    #snapshots = ["benchM", "NsPNG_EDE_F500","NsPNG_EDE_F1833", "NG_ViVi","NG_Fminus500_ViVi"]
    #labels = [r"$\Lambda$CDM", r"$f_{\rm NL} = -300~\&~{\rm EDE}$",  r"$f_{\rm NL} = -1100~\&~{\rm EDE}$", r"$f_{\rm NL} = -500~\&~{\rm mixed~DM}$", r"$f_{\rm NL} = 500~\&~ {\rm mixed~DM}$"]


    #lss = ["-", ":", "--",  "--","--"]
    #couleurs = ["blue", "darkred", "darkred","darkorange","violet"]


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


    outer = gridspec.GridSpec(nrows=1, ncols=1)

    axs = []
    for row in range(1):
        for col in range(1):
            inner = gridspec.GridSpecFromSubplotSpec(nrows=2, ncols=1, subplot_spec=outer[row, col], hspace=0)
            axs += [plt.subplot(cell) for cell in inner]


    for i in range(len(snapshots)):
        k = 1
        moyennes = []
        err = []

        ls = lss[i]
        couleur = couleurs[i]
        label = labels[i]
        moyennes_lcdm =[]
        err_lcdm = []
    
        for j in indices_z:

            long_lcdm = np.load(f"/data100/fcastillo/RESULT/{snapshots[0]}/{k}_densite_smooth2_c0.1_len_fil.txt.npy")
            try : long = np.load(f"/data100/fcastillo/RESULT/{snapshots[i]}/{k}_densite_smooth2_c0.1_len_fil.txt.npy")
            except :     long = np.load(f"/data100/fcastillo/RESULT/{snapshots[i]}/{j}_densite_smooth2_c0.1_len_fil.txt.npy")
            moyennes.append(np.mean(long))
            err.append(1/sqrt(len(long)) *np.std(long))
            moyennes_lcdm.append(np.mean(long_lcdm))
            err_lcdm.append(1/sqrt(len(long_lcdm)) *np.std(long_lcdm))
            k += 1
        

        axes = plt.gca()

        a = 1/(1+np.array([3,1,0.25,0]))

        
        moyennes_lcdm = np.array(moyennes_lcdm)
        err_lcdm = np.array(err_lcdm)
        err= np.array(err)

        for d in range(2):

            moyennes = np.array(moyennes)

            #plt.scatter(np.log(np.array([32,3,1,0.25,0][:len(moyennes)])), moyennes, color=couleur)


            try : err_ratio = np.sqrt((moyennes_lcdm**2*err**2 + moyennes**2*err_lcdm**2)/moyennes_lcdm**4)
            except : 
                print(err,err_lcdm)
                err_ratio = np.array([0,0,0,0])

            if d == 0 : axs[d].errorbar(np.array([3,1,0.25,0]), moyennes/500**3,ls=ls, color=couleur, label=label, yerr = err/500**3)
            if d == 1 : axs[d].errorbar(np.array([3,1,0.25,0]), (moyennes-moyennes_lcdm)/moyennes_lcdm,ls=ls, color=couleur, label=label, yerr = err_ratio)

            #print(snapshots)
            if d == 0 : 
                #print(snapshots[i]+" & ",end="")
                for z in [0,1,3]:
                    print(str(round(100*((moyennes-moyennes_lcdm)/moyennes_lcdm)[z],3))+" \%",end="")
                    print(" & ",end="")
                    #else : print(" \\\\")

            axes.set_xlabel(r"$z$")
            if d == 0 : axs[d].set_ylabel(r"${\rm Mean~length}~/~V~~[{\rm Mpc / h}]^{-2}$")
            if d == 1 : axs[d].set_ylabel(r"$\Delta / \Lambda$CDM")

            if d == 0 : axs[d].legend(fontsize=8)

    axs[0].invert_xaxis()
    axs[1].invert_xaxis()
    plt.savefig(f"len_moyenne_EDE.pdf")
    plt.savefig(f"len_moyenne_EDE.png")

