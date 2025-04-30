import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from astropy.io import fits
from matplotlib import gridspec



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

Redshifts = [3,1,0.25,0]
pre = "/data100/fcastillo/RESULT/"



for i in range(4):

    outer = gridspec.GridSpec(nrows=3, ncols=3)

    plt.figure(figsize=(9,9))


    axs = []
    for row in range(3):
        for col in range(3):
            inner = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=1, subplot_spec=outer[row, col], hspace=0)
            axs += [plt.subplot(cell) for cell in inner]


    k = 0

    data = pre + snapshots[0]+"/"+str(i)+"_densite_smooth2.fits"

    hdul = fits.open(data)
    lcdm = hdul[0].data
    hdul.close()

    sum_ = np.sum(lcdm[:,:,0:10],axis=2)
    mean_ = np.mean(sum_)
    std_ = np.std(sum_)
    z0 = mean_ - 2*std_     
    z1 = mean_ +2*std_


    plt.title(r"$z = "+str(Redshifts[i])+r"$")
    for n in indices_hdm: 
        k +=1 
        if n <= 8 : redshifts = range(1,5)
        else : redshifts = indices_z
        z = redshifts[i]
        print(z)
        axes = axs[k-1]
        axes.title.set_text(labels[n])


        data = pre + snapshots[n]+"/"+str(z)+"_densite_smooth2.fits"

        hdul = fits.open(data)
        field = hdul[0].data
        hdul.close()


        if n < 8:
            im = axes.imshow(np.sum(field[:,:,0:10],axis=2), origin="lower",vmin = z0, vmax = z1)
        else :
            im = axes.imshow(np.sum(field[:,:,0:10],axis=0), origin="lower",vmin = z0, vmax = z1)

        if (k) % 3 == 0 : plt.colorbar(im, ax=axes)
        if not k-1 > 5 : axes.xaxis.set_visible(False)  
        if not (k-1) % 3 == 0 : axes.yaxis.set_visible(False)  

        axes.set_xlim(0,500)
        axes.set_ylim(0,500)

        axes.set_xlabel(r"$\rm X [Mpc / h]$")
        axes.set_ylabel(r"$\rm Y [Mpc / h]$")

    plt.tight_layout()
    plt.savefig(f"field_{Redshifts[i]}.pdf")
    plt.clf()