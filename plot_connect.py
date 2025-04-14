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

    #snapshots = ["benchM", "G_ViVi","NG_F500", "NG_Fminus500","NEDE"]
    #labels = [r"$\Lambda$CDM", r"$m_{\rm WDM} = 10  {\rm eV}, f_{\rm WDM} = 2 \%$",  r"$f_{\rm NL} = -500$", r"$f_{\rm NL} = 500$ ", r"${\rm EDE}$"]


    #lss = ["-", "-.", "-",  "-","-"]
    #couleurs = ["blue", "green", "darkorange","violet","darkred"]



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

    z= [15,12, 10, 8, 5,3,1,0.5,0.25,0]
    indices_z = [5,6,8,9]
    


    outer = gridspec.GridSpec(nrows=2, ncols=2)

    axs = []
    for row in range(2):
        for col in range(2):
            inner = gridspec.GridSpecFromSubplotSpec(nrows=2, ncols=1, subplot_spec=outer[row, col], hspace=0)
            axs += [plt.subplot(cell) for cell in inner]


    k = 0
    for i in [1,2,3,4]:
        k += 1
        z_k = indices_z[k-1]
        lcdm = np.load(f"/data100/fcastillo/RESULT/{snapshots[0]}/{i}_densite_smooth2_c0.1_connect_fil.txt.npy")
        hist_lcdm = np.histogram(lcdm,  density= True, range = [0, 10], bins=nbins)
        hist_lcdm = hist_lcdm[0]


        for d in range(2):
            place = places[str(i) + str(d)]

            print(axs)

            axes = axs[place-1]

            if d == 0 : 
                axes.title.set_text (f"z = {Redshifts[i]}")


            for j in range(len(snapshots)):

                ls = lss[j]
                couleur = couleurs[j]
                label = labels[j]


                try :
                        
                    try : connect = np.load(f"/data100/fcastillo/RESULT/{snapshots[j]}/{i}_densite_smooth2_c0.1_connect_fil.txt.npy")
                    except : connect = np.load(f"/data100/fcastillo/RESULT/{snapshots[j]}/{z_k}_densite_smooth2_c0.1_connect_fil.txt.npy")
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

                    if j == len(snapshots)-1 and i == 1 and d == 0: 
                        axes.legend() 
                
                except : pass

    if i == 0:
        plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)

    plt.savefig(f"connect_{nbins}_EDE.pdf")
    plt.savefig(f"connect_{nbins}_EDE.png")

    plt.clf()

    plt.figure()


    outer = gridspec.GridSpec(nrows=1, ncols=1)

    axs = []
    for row in range(1):
        for col in range(1):
            inner = gridspec.GridSpecFromSubplotSpec(nrows=2, ncols=1, subplot_spec=outer[row, col], hspace=0)
            axs += [plt.subplot(cell) for cell in inner]

    for i in range(len(snapshots)):
        moyennes = []
        err = []

        ls = lss[i]
        couleur = couleurs[i]
        label = labels[i]
        moyennes_lcdm = []
        err_lcdm = []
     
        for j in range(1,5):
            z_k = indices_z[j-1]
                    
            connect_lcdm = np.load(f"/data100/fcastillo/RESULT/{snapshots[0]}/{j}_densite_smooth2_c0.1_connect_fil.txt.npy")
            try : connect = np.load(f"/data100/fcastillo/RESULT/{snapshots[i]}/{j}_densite_smooth2_c0.1_connect_fil.txt.npy")
            except : connect = np.load(f"/data100/fcastillo/RESULT/{snapshots[i]}/{z_k}_densite_smooth2_c0.1_connect_fil.txt.npy")
            moyennes.append(np.mean(connect))
            moyennes_lcdm.append(np.mean(connect_lcdm))
            err.append(1/sqrt(len(connect)) *np.std(connect))
            err_lcdm.append(1/sqrt(len(connect_lcdm)) *np.std(connect_lcdm))
            
            

        print(moyennes)
        axes = plt.gca()

        if len(moyennes )== 5:
            a = 1/(1+np.array([32,3,1,0.25,0]))
        else :
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

            if d == 0 : axs[d].errorbar(np.array([3,1,0.25,0]), moyennes,ls=ls, color=couleur, label=label, yerr = err)
            if d == 1 : axs[d].errorbar(np.array([3,1,0.25,0]), (moyennes-moyennes_lcdm)/moyennes_lcdm,ls=ls, color=couleur, label=label, yerr = err_ratio)

            axes.set_xlabel(r"$z$")
            if d == 0 : axs[d].set_ylabel("Mean connectivity")
            if d == 1 : axs[d].set_ylabel(r"$\Delta / \Lambda$CDM")

            if d == 0 : axs[d].legend(fontsize=8)

    axs[0].invert_xaxis()
    axs[1].invert_xaxis()

    plt.savefig(f"connect_moyenne_EDE.pdf")
    plt.savefig(f"connect_moyenne_EDE.png")

