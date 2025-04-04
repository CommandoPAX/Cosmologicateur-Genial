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
   
    snapshots = ["benchM","NG_F500","G_m500","NG_F500_m500","NG_Fminus500","NG_Fminus500_m500", "G_ViVi", "NG_ViVi", "NG_Fminus500_ViVi"]
    labels = ["LCDM", "fnl = -500", "m = 500 eV", "WDM & fnl = -500", "fnl = 500", "WDM & fnl = 500", r"$m_{\rm WDM} = 10$ ev, $f_{\rm WDM}$ = 2%", "fnl = -500 & mixed DM", "fnl = 500 & mixed DM"]

    indices_hdm = [0,1,4,6,7,8]
    indices_hdm = [0,2,6] #WDM !

    ls = ["-", "-", "-.", "--", "-", "--", "-", "--", "--"]
    couleurs = ["blue", "orange", "green", "orange", "fuchsia", "fuchsia", "green", "orange", "fuchsia"]

    fig = plt.figure(figsize=(12,6))
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

    R = 2
    P = 5
    

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



            for d in range(2):
                place = places[str(p) + str(d)]

                print(axs)

                axes = axs[(place-1)+(i//2-1)*8]

                if d == 0 : 
                    axes.title.set_text (rf"{["P","F","W","V"][p]}{["P","F","W","V"][p]},  $z = $"+str(Redshifts[i]))


                for j in indices_hdm:
                    if True:
                        lcdm = np.load(f"/data100/fcastillo/RESULT/extrema/snapshot_{0}_{i}_zeta_{p}_{p}_s{R}_P{P}.txt.npy")

                        zeta = np.load(f"/data100/fcastillo/RESULT/extrema/snapshot_{j}_{i}_zeta_{p}_{p}_s{R}_P{P}.txt.npy")
                        #zeta[0] = 0

                        r_small = np.linspace(R, 10, 80)  # 10 points entre 0 et 1
                        r_large = np.geomspace(10, 40, 20)  # 30 points entre 1 et 40 (logarithmique)
                        r_bins = np.concatenate((r_small, r_large))


                        #print(np.shape(data[p]))
                        #print(data)
                        if d == 0 : 
                            axes.plot(r_bins[1:], zeta, color=couleurs[j], ls=ls[j],label=labels[j])
                            if p == 0 :axes.set_ylabel(r"$1 + \zeta (r)$")
                            axes.xaxis.set_visible(False)
                        else : 
                            axes.plot(r_bins[1:], zeta - lcdm, color=couleurs[j], ls=ls[j],label=labels[j])
                            if p == 0 : axes.set_ylabel(r"$\Delta$")
                        if d ==1 and i == 4: axes.set_xlabel("r [Mpc / h]")
                        if j == 6 and i == 2 and d == 0 and p == 0: 
                            axes.legend(fontsize=8) 
                        axes.set_xlim(0,10)

        if i == 0:
            plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)

        plt.tight_layout()
        plt.savefig(f"corr_auto.pdf")
        plt.savefig(f"corr_auto.png")

        #plt.clf()

    plt.clf()
    fig = plt.figure(figsize=(8,6))
    #plt.title(  rf"v$_{p}$")
    plt.axis("off")
    #plt.tight_layout()

    outer = gridspec.GridSpec(nrows=2, ncols=2)

    axs = []
    for row in range(2):
        for col in range(2):
            inner = gridspec.GridSpecFromSubplotSpec(nrows=2, ncols=1, subplot_spec=outer[row, col], hspace=0)
            axs += [plt.subplot(cell) for cell in inner]

    for i in [2,4] :



        for p in range(2):


            _type = ["PF", "VW"] [p]

            for d in range(2):
                place = places[str(p) + str(d)]

                print(axs)

                axes = axs[(place-1)+(i//2 -1)*4]

                if d == 0 : 
                    axes.title.set_text (rf"{_type},  $z = $"+str(Redshifts[i]))


                for j in indices_hdm:
                    if True:
                        if p == 0 : 
                            a = 1
                            b = 0
                        elif p == 1 :
                            a = 3
                            b = 2
                        lcdm = np.load(f"/data100/fcastillo/RESULT/extrema/snapshot_{0}_{i}_zeta_{a}_{b}_s{R}_P{P}.txt.npy")

                        zeta = np.load(f"/data100/fcastillo/RESULT/extrema/snapshot_{j}_{i}_zeta_{a}_{b}_s{R}_P{P}.txt.npy")
                        #zeta[0] = 0

                        r_small = np.linspace(R, 10, 80)  # 10 points entre 0 et 1
                        r_large = np.geomspace(10, 40, 20)  # 30 points entre 1 et 40 (logarithmique)
                        r_bins = np.concatenate((r_small, r_large))


                        #print(np.shape(data[p]))
                        #print(data)
                        if d == 0 : 
                            axes.plot(r_bins[1:], zeta, color=couleurs[j], ls=ls[j],label=labels[j])
                            axes.set_ylabel(r"$1 + \zeta (r)$")
                            axes.xaxis.set_visible(False)
                            axes.set_ylim(-0.2,1)
                            if p == 1 : axes.set_ylim(0,3)
                        else : 
                            axes.plot(r_bins[1:], zeta - lcdm, color=couleurs[j], ls=ls[j],label=labels[j])
                            axes.set_ylabel(r"$\Delta$")
                        if d ==1 and i == 4: axes.set_xlabel("r [Mpc / h]")
                        if j == 6 and i == 2 and d == 0 and p == 0: 
                            axes.legend(loc="upper left",fontsize=8) 

                        axes.set_xlim(0,10)


        if i == 0:
            plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)

        plt.tight_layout()
        plt.savefig(f"corr_PF_VW.pdf")
        plt.savefig(f"corr_PF_VW.png")

        #plt.clf()

    plt.clf()
    fig = plt.figure(figsize=(12,6))
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


            _type = ["PW", "PV", "FW","FV"] [p]

            for d in range(2):
                place = places[str(p) + str(d)]

                print(axs)

                axes = axs[(place-1)+(i//2 -1)*8]

                if d == 0 : 
                    axes.title.set_text (rf"{_type},  $z = $"+str(Redshifts[i]))


                for j in indices_hdm:
                    if True:
                        if p == 0 : 
                            a = 2
                            b = 0
                        elif p == 1 :
                            a = 3
                            b = 0
                        elif p == 2 :
                            a = 2
                            b = 1
                        elif p == 3 :
                            a = 3
                            b = 1
                        lcdm = np.load(f"/data100/fcastillo/RESULT/extrema/snapshot_{0}_{i}_zeta_{a}_{b}_s{R}_P{P}.txt.npy")

                        zeta = np.load(f"/data100/fcastillo/RESULT/extrema/snapshot_{j}_{i}_zeta_{a}_{b}_s{R}_P{P}.txt.npy")
                        #zeta[0] = 0

                        r_small = np.linspace(R, 10, 80)  # 10 points entre 0 et 1
                        r_large = np.geomspace(10, 40, 20)  # 30 points entre 1 et 40 (logarithmique)
                        r_bins = np.concatenate((r_small, r_large))


                        #print(np.shape(data[p]))
                        #print(data)
                        if d == 0 : 
                            axes.plot(r_bins[1:], zeta, color=couleurs[j], ls=ls[j],label=labels[j])
                            if p == 0 : axes.set_ylabel(r"$1 + \zeta (r)$")
                            axes.xaxis.set_visible(False)
                            #axes.set_ylim(-1,1)
                        else : 
                            axes.plot(r_bins[1:], zeta - lcdm, color=couleurs[j], ls=ls[j],label=labels[j])
                            if p == 0 :axes.set_ylabel(r"$\Delta$")
                        if d ==1 and i == 4: axes.set_xlabel("r [Mpc / h]")
                        if j == 6 and i == 2 and d == 0 and p == 0: 
                            axes.legend(fontsize=8) 
                        axes.set_xlim(0,20)

        #if i == 0:
        #    plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)

        plt.tight_layout()
        plt.savefig(f"corr_autres.pdf")
        plt.savefig(f"corr_autres.png")

        #plt.clf()
