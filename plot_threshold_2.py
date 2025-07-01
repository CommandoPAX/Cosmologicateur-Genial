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
labels = [r"$\Lambda$CDM", "fnl = -500", "m = 500 eV", "WDM & fnl = -500", "fnl = 500", "WDM & fnl = 500",r"$m_{\rm WDM} = 10$ ev, $f_{\rm WDM}$ = 2%"]

matplotlib.rcParams.update({'font.size': 18})
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams["figure.facecolor"]='w'
matplotlib.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'

indices_z = [5,6,9]

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

    snapshots = ["benchM","NG_F500","G_m500","NG_F500_m500","NG_Fminus500","NG_Fminus500_m500", "G_ViVi","NG_ViVi" , "NG_Fminus500_ViVi","NEDE", "NsPNG_EDE_F500","NBM"]

    labels = [r"$\Lambda{\rm CDM}$", "fnl = -500", "m = 500 eV", "WDM & fnl = -500", "fnl = 500", "WDM & fnl = 500", r"mixed DM", r"$f_{\rm NL}^0 = -500$ \& mixed DM", r"$f_{\rm NL}^0$ = 500 \& mixed DM", r"EDE", r"$f_{\rm NL}^0 = -300$ \& EDE", r"$f_{\rm NL}^0 = -1100$ \& EDE",r"$\Lambda {\rm CDM~2}$"]

    indices_hdm = [0,1,4,6,7,8]
    indices_hdm = [0,1,4,6,9]
    indices_hdm = [0,12]

    #indices_hdm = [0,2,6] #WDM !

    Points = ["P","F","W","V"]


    lss = ["-", "-", "-.", "--", "-", "--", "-", "--", "--"]
    couleurs = ["blue", "orange", "green", "orange", "fuchsia", "fuchsia", "green", "orange", "fuchsia"]

    snapshots = ["benchM","NG_F500","G_m500","NG_F500_m500","NG_Fminus500","NG_Fminus500_m500", "G_ViVi","NG_ViVi" , "NG_Fminus500_ViVi","NEDE", "NsPNG_EDE_F500","NBM"]


    labels = [r"$\Lambda{\rm CDM}$", r"$f_{\rm NL}^0 = -500$", "m = 500 eV", "WDM & fnl = -500", r"$f_{\rm NL}^0 = 500$", "WDM & fnl = 500", r"${\rm mixed~DM}$", r"$f_{\rm NL}^0 = -500~\&~{\rm mixed~DM}$", r"$f_{\rm NL}^0 = 500~\&~{\rm mixed~DM}$", r"${\rm EDE}$", r"$f_{\rm NL}^0 = -300~\&~{\rm EDE}$", r"$f_{\rm NL}^0 = -1100~\&~ {\rm EDE}$",r"$\Lambda {\rm CDM~2}$"]

    #indices_hdm = [0,7,8,10,11]
    #indices_hdm = [0,2,6] #WDM !

    lss = ["-", "-", "-.", "--", "-", "--", "-.", "--", "--","-",":","--","--"]
    couleurs = ["blue", "darkorange", "green", "darkorange", "violet", "violet", "green", "darkorange", "violet","darkred","darkred","darkred","blue"]



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

    outer = gridspec.GridSpec(nrows=3, ncols=4)
    R = 2

    axs = []
    for row in range(3):
        for col in range(4):
            inner = gridspec.GridSpecFromSubplotSpec(nrows=2, ncols=1, subplot_spec=outer[row, col], hspace=0)
            axs += [plt.subplot(cell) for cell in inner]

    k = 0
    for i in [1,2,4] :
        z_k = indices_z[k]
        k += 1


        for p in range(4):
            lcdm = np.load(f"/data100/fcastillo/RESULT/extrema/snapshot_{0}_{i}_threshold_s{R}.txt.npy")



            for d in range(2):
                place = places[str(p) + str(d)]

                print(axs)

                axes = axs[(place-1)+(min(i-1,2))*8]

                if d == 0 : 
                    axes.title.set_text (r"$\mathcal{"+rf"{Points[p]}"+r"},  z = "+str(Redshifts[i])+r"$")


                for j in indices_hdm:
                    try : count = np.load(f"/data100/fcastillo/RESULT/extrema/snapshot_{j}_{i}_threshold_s{R}.txt.npy")
                    except: count = np.load(f"/data100/fcastillo/RESULT/extrema/snapshot_{j}_{z_k}_threshold_s{R}.txt.npy")
                    #zeta[0] = 0
                    Npoints = len(count)
                    print(Npoints)
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
                    if i == 2 and p == 3 : axes.set_xlim(-6,0)
                    if i == 4 and p == 3 : axes.set_xlim(-6,0)
                    if i == 2 and p == 2 : axes.set_xlim(-3,2)
                    if i == 4 and p == 2 : axes.set_xlim(-3,2)

                    if d == 0 : 
                        axes.plot(X, count[:,p] , color=couleur, ls=ls,label=label)
                        axes.xaxis.set_visible(False)
                        if p == 0 : axes.set_ylabel(r"$\frac{1}{N} \frac{dN}{d\nu}$")

                    else : 
                        mask = lcdm[:,p] > 0.1 * np.std(lcdm[:,p])
                        axes.plot(X, ((count[:,p]-lcdm[:,p])/np.max(lcdm[:,p])), color=couleur, ls=ls,label=label)
                        if p == 0 : axes.set_ylabel(r"$\Delta / {\rm max}_{\Lambda {\rm CDM}}$")
                    if d ==1 and i == 4: 
                        axes.set_xlabel(r"$\nu [\sigma]$")
                        #axes.set_xlim(-6,6)

                    if j == indices_hdm[len(indices_hdm)-1] and i == 1 and d == 0 and p == 0: 
                        pass#axes.legend(fontsize = 8) 
                #except: pass

        if i == 0:
            plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)


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
        fig.legend(unique_handles, unique_labels, loc='upper center', ncol=3, frameon=True,fontsize=18)

    plt.tight_layout(rect=[0, 0, 1, 0.9])
    plt.savefig(f"crit_threshold_s{R}.pdf")
    plt.savefig(f"crit_threshold_s{R}.png")
