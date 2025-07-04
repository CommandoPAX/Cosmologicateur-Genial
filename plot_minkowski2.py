import numpy as np
import os
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from astropy.io import fits
from math import*
from matplotlib import gridspec
import matplotlib

snapshots = ["benchM","NG_F500","G_m500","NG_F500_m500","NG_Fminus500","NG_Fminus500_m500","G_ViVi"]
labels = [r"$\Lambda{\rm CDM}$", "fnl = -500", "m = 500 eV", "WDM & fnl = -500", "fnl = 500", "WDM & fnl = 500",r"$m_{\rm WDM} = 10$ ev, $f_{\rm WDM}$ = 2%"]

z= [15,12, 10, 8, 5,3,1,0.5,0.25,0]
indices_z = [5,6,9]

matplotlib.rcParams.update({'font.size': 18})
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

    #snapshots = ["benchM","NG_F500","G_m500","NG_F500_m500","NG_Fminus500","NG_Fminus500_m500","G_ViVi"]
    #labels = [r"$\Lambda$CDM", "fnl = -500", "m = 500 eV", "WDM & fnl = -500", "fnl = 500", "WDM & fnl = 500",r"$m_{\rm WDM} = 10$ ev, $f_{\rm WDM}$ = 2%"]

    #snapshots = ["benchM","NG_F500","G_ViVi","NG_ViVi","NG_Fminus500","NG_Fminus500_ViVi"]
    #labels = [r"$\Lambda$CDM", "fnl = -500", r"$m_{\rm WDM} = 10$ ev, $f_{\rm WDM}$ = 2%", "fnl = -500 & mixed DM", "fnl = 500", "fnl = 500 & mixed DM"]


    #ls = ["-", "-", "-.", "--", "-", "--","-"]
    #couleurs = ["blue", "orange", "green", "orange", "fuchsia", "fuchsia"]


    snapshots = ["benchM", "NsPNG_EDE_F500","NsPNG_EDE_F1833", "NG_ViVi","NG_Fminus500_ViVi"]
    labels = [r"$\Lambda {\rm CDM}$", r"$f_{\rm NL}^0 = -300~\&~ \rm{EDE}$",  r"$f_{\rm NL}^0 = -1100~\&~{\rm EDE}$", r"$f_{\rm NL}^0 = -500~\&~{\rm mixed~DM}$", r"$f_{\rm NL}^0 = 500~\&~{\rm mixed~DM}$"]


    ls = ["-", ":", "--",  "--","--"]
    couleurs = ["blue", "darkred", "darkred","darkorange","violet"]

    snapshots = ["benchM", "G_ViVi","NG_F500", "NG_Fminus500","NEDE"]
    labels = [r"$\Lambda{\rm CDM}$", r"${\rm mixed~DM}$",  r"$f_{\rm NL}^0 = -500$", r"$f_{\rm NL}^0 = 500$ ", r"${\rm EDE}$"]


    ls = ["-", "-.", "-",  "-","-"]
    couleurs = ["blue", "green", "darkorange","violet","darkred"]

    snapshots = ["benchM", "G_m500","G_ViVi"]
    labels = [r"$\Lambda{\rm CDM}$", r"$m_{\rm WDM} = 500 {\rm ~eV}$", r"${\rm mixed~DM}$"]


    ls = ["-", "-", "-."]
    couleurs = ["blue", "green", "green"]

    snapshots = ["benchM", "NBM"]
    labels = [r"$\Lambda{\rm CDM}$", r"$\Lambda{\rm CDM~2}$"]


    ls = ["-", "--"]
    couleurs = ["blue", "blue"]


    plt.figure(figsize=(15,10))
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

    outer = gridspec.GridSpec(nrows=3, ncols=4)

    axs = []
    for row in range(3):
        for col in range(4):
            inner = gridspec.GridSpecFromSubplotSpec(nrows=2, ncols=1, subplot_spec=outer[row, col], hspace=0)
            axs += [plt.subplot(cell) for cell in inner]
    k = 0

    for i in [1,2,4] :

        z_k = indices_z[k]
        k +=1

        for p in range(4):
            lcdm = np.load(f"/home/fcastillo/minkowski/minkowski_{0}_{i}.txt.npy")
            if i in [2,4]:lcdmzoom = np.load(f"/data100/fcastillo/RESULT/benchM/{i}_minkowski_zoom.txt.npy")

            X = np.linspace(-3,3,61)
            if i in [2,4]: 
                X1 = X[X<-1]
                X2 = np.linspace(-1,1,101)
                X3 = X[X>1]

                X = np.concatenate((X1,X2,X3))

            if p == 0 and i in [2,4] : X = X2

            for d in range(2):
                place = places[str(p) + str(d)]

                axes = axs[(place-1)+(min(i-1,2))*8]

                if d == 0 : 
                    axes.title.set_text (rf"$v_{p}$,  $z = "+str(Redshifts[i])+"$")


                for j in range(len(snapshots)):
                    if True:
                        lcdm = np.load(f"/home/fcastillo/minkowski/minkowski_{0}_{i}.txt.npy")

                        try:
                            data = np.load(f"/home/fcastillo/minkowski/minkowski_{j}_{i}.txt.npy")
                        except :
                            try : data= np.load(f"/data100/fcastillo/RESULT/{snapshots[j]}/{i}_minkowski.txt.npy")
                            except: data= np.load(f"/data100/fcastillo/RESULT/{snapshots[j]}/{z_k}_minkowski.txt.npy")
                        if i in [2,4]:
                            try :
                                datazoom = np.load(f"/data100/fcastillo/RESULT/{snapshots[j]}/{i}_minkowski_zoom.txt.npy")
                            except :
                                datazoom = np.load(f"/data100/fcastillo/RESULT/{snapshots[j]}/{z_k}_minkowski_zoom.txt.npy")
                            data1 = data[:,:20]
                            data2 = data[:,41:]

                            data = np.concatenate([data1,datazoom, data2],axis=1)


                            lcdm = np.concatenate([lcdm[:,:20],lcdmzoom, lcdm[:,41:]],axis=1)

                            if p == 0 :
                                data = datazoom
                                lcdm = lcdmzoom



                        if p == 0 and i in [2,4] : mask = (X > -1) & (X <1)
                        elif p!=0 and i in [2,4] : mask = (X > -1)&(X < 3)
                        else : mask = (X>-3)&(X<3)

                        if d == 1 : 
                            max_delta = np.argmax(np.abs(data[p][mask]-lcdm[p][mask]))
                            if j !=0 : print(labels[j]+f", v{p}, z = "+str(Redshifts[i])+r" : $\sim$"+str(round(100*((data[p][mask][max_delta]-lcdm[p][mask][max_delta])/np.max(np.abs(lcdm[p][mask]))),1))+" \%")
                            #else : print(" \\\\")

                        if d == 0 : 
                            axes.plot(X[mask], data[p][mask], color=couleurs[j], ls=ls[j],label=labels[j])
                            axes.set_ylabel(rf"$v_{p}$")
                            axes.xaxis.set_visible(False)
                        else :
                            axes.plot(X[mask], ((data[p][mask] - lcdm[p][mask])), color=couleurs[j], ls=ls[j],label=labels[j])
                            axes.set_ylabel(r"$\Delta$")
                        if d ==1 and i == 4: axes.set_xlabel(r"${\rm threshold}~[\sigma]$")
                        if j == len(snapshots)-1 and i == 1 and d == 0 and p == 0: 
                            pass#axes.legend(fontsize=12) 

                        if p == 0 and i in [2,4] : axes.set_xlim(-1,1)
                        elif p!=0 and i in [2,4] : axes.set_xlim(-1,3)
                        else : axes.set_xlim (-3,3)

        if i == 0:
            plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)

        # Récupérer tous les handles/labels
        handles = []
        labels_all = []
        ax = axs[0]
        h, l = ax.get_legend_handles_labels()
        handles.extend(h)
        labels_all.extend(l)

        # Éliminer les doublons (en conservant l'ordre)
        seen = set()
        unique_handles_labels = [(h, l) for h, l in zip(handles, labels_all) if not (l in seen or seen.add(l))]
        unique_handles, unique_labels = zip(*unique_handles_labels)

        # Légende globale
        fig = plt.gcf()
        fig.legend(unique_handles, unique_labels, loc='upper center', ncol=5, frameon=True,fontsize=18)

        # Ajuster la mise en page pour laisser la place à la légende


        plt.tight_layout(rect=[0, 0, 1, 0.95])
        plt.savefig(f"v_tout.pdf")
        plt.savefig(f"v_tout.png")

        #plt.clf()

