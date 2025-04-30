import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from astropy.io import fits


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

redshifts = [3,1,0.25,0]

k = 0
for i in range(4):
    plt.figure()
    plt.title(r"z = "+str(redshifts[i]))
    for n in indices_hdm: 
        k +=1 
        if n <= 8 : redshifts = range(1,5)
        else : redshifts = indices_z
        for i in redshifts:
            z = redshifts[i]
            plt.subplot(int("33"+str(k)))
            axes = plt.gca()
            axes.title.set_text(labels[n])

            pre = "/data100/fcastillo/RESULT/"

            data = pre + snapshots[n]+"/"+str(i)+"_densite_smooth2.fits"

            hdul = fits.open(data)
            field = hdul[0].data
            hdul.close()


            im = axes.imshow(np.sum(field,axis=2), origin="lower", vmax = 10, vmin = -2)
            axes.colorbar(im)

            axes.set_xlim(0,10)
            axes.set_ylim(0,10)

            axes.set_xlabel(r"$\rm X [Mpc / h]$")
            axes.set_ylabel(r"$\rm Y [Mpc / h]$")

    plt.savefig(f"field_{redshifts[i]}.pdf")