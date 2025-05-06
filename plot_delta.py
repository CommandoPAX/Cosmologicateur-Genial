import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from astropy.io import fits
from matplotlib import gridspec
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable


matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams["figure.facecolor"]='w'
matplotlib.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'

snapshots = ["benchM","NG_F500","G_m500","NG_F500_m500","NG_Fminus500","NG_Fminus500_m500", "G_ViVi","NG_ViVi" , "NG_Fminus500_ViVi","NEDE", "NsPNG_EDE_F500","NsPNG_EDE_F1833"]
labels = [r"$\Lambda$CDM", r"$f_{\rm NL} = -500$", "m = 500 eV", "WDM & fnl = -500", r"$f_{\rm NL} = 500$", "WDM & fnl = 500", r"$m_{\rm WDM} = 10~{\rm eV}, f_{\rm WDM} = 2~\%$", r"$f_{\rm NL} = -500~\&~{\rm mixed~DM}$", r"$f_{\rm NL} = 500~\&~{\rm mixed~DM}$", r"${\rm EDE}$", r"$f_{\rm NL} = -300~\&~{\rm EDE}$", r"$f_{\rm NL} = -1100~\&~{\rm EDE}$"]

z= [15,12, 10, 8, 5,3,1,0.5,0.25,0]
indices_z = [5,6,8,9]

indices_hdm = [0,1,4,6,7,8,9,10,11]
#indices_hdm = [0,2,6]
#indices_hdm = [0,1,4,6,9]
indices_hdm = [0,1,4,6,9]

ls = ["-", "-", "-.", "--", "-", "--", "-.", "--", "--","-",":","--"]
couleurs = ["blue", "darkorange", "green", "darkorange", "violet", "violet", "green", "darkorange", "violet","darkred","darkred","darkred"]

Redshifts = [3,1,0.25,0]
pre = "/data100/fcastillo/RESULT/"



if True : #for Xmax in [50,100,500]:
    for i in [0] :#range(4):

        outer = gridspec.GridSpec(nrows=3, ncols=3)

        plt.figure(figsize=(6,6))


        k = 0

        data = pre + snapshots[0]+"/"+str(i+1)+"_densite_smooth2.fits"


        plt.title(r"$z = "+str(Redshifts[i])+r"$")

        plt.subplot(211)
        for n in indices_hdm: 
            k +=1 
            if n <= 8 : redshifts = range(1,5)
            else : redshifts = indices_z
            z = redshifts[i]
            print(z)
            axes = plt.gca()


            data = pre + snapshots[n]+"/"+str(z)+"_densite_smooth2.fits"

            hdul = fits.open(data)
            field = hdul[0].data
            hdul.close()
            field = field.reshape(-1, 1) 
            hist = np.histogram(field,  density= True,bins=1000,range = [-2,3])
            hist = hist[0]

            plt.plot(np.arange(-2,3,5 / 1000),hist,color=couleurs[n], ls=ls[n],label=labels[n])

        plt.subplot(212)

        data = pre + snapshots[0]+"/"+str(range(1,5)[i])+"_densite_smooth2.fits"

        hdul = fits.open(data)
        lcdm = hdul[0].data
        hdul.close()
        lcdm = lcdm.reshape(-1,1)

        hist_lcdm = np.histogram(lcdm,  density= True,bins=1000,range = [-2,3])
        hist_lcdm = hist[0]


        for n in indices_hdm: 
            k +=1 
            if n <= 8 : redshifts = range(1,5)
            else : redshifts = indices_z
            z = redshifts[i]
            print(z)
            axes = plt.gca()


            data = pre + snapshots[n]+"/"+str(z)+"_densite_smooth2.fits"

            hdul = fits.open(data)
            field = hdul[0].data
            hdul.close()
            field = field.reshape(-1, 1) 
            hist = np.histogram(field,  density= True,bins=1000,range = [-2,3])
            hist = hist[0]

            plt.plot(np.arange(-2,3,5 / 1000),hist-hist_lcdm,color=couleurs[n], ls=ls[n],label=labels[n])


        plt.legend()
        plt.tight_layout()
        plt.savefig(f"delta.pdf")
        plt.clf()
