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

indices_hdm = [0,1,4,6,7,8]
#indices_hdm = [0,2,6]
#indices_hdm = [0,1,4,6,9]

Redshifts = [3,1,0.25,0]
pre = "/data100/fcastillo/RESULT/"



for minmax in range(2) :
    ims = []

    for Xmax in [50]:
        for i in range(4):

            outer = gridspec.GridSpec(nrows=2, ncols=3)

            plt.figure(figsize=(9,6))


            axs = []
            for row in range(2):
                for col in range(3):
                    inner = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=1, subplot_spec=outer[row, col], hspace=0)
                    axs += [plt.subplot(cell) for cell in inner]


            k = 0

            data = pre + snapshots[0]+"/"+str(i+1)+"_densite_smooth2.fits"

            hdul = fits.open(data)
            lcdm = hdul[0].data
            hdul.close()


            if minmax == 0 : max_index_flat = np.argmax(lcdm)
            else : max_index_flat = np.argmin(lcdm)
            max_position = np.unravel_index(max_index_flat, lcdm.shape)

            X0 = max_position[0]
            Y0 = max_position[1]
            Z0 = max_position[2]

            sum_lcdm = np.sum(lcdm[X0-25:X0+25,Y0-25:Y0+25,Z0-1:Z0+1],axis=2)
            rho_m = 1073741824000000
            #mass = rho_m * sum_ + rho_m
            


            plt.title(r"$z = "+str(Redshifts[i])+r"$")
            for n in indices_hdm: 
                k +=1 
                if n <= 8 : redshifts = range(1,5)
                else : redshifts = indices_z
                z = redshifts[i]
                print(z)
                axes = axs[k-1]
                axes.title.set_text(labels[n]+r"$ - \Lambda{\rm CDM}$")


                data = pre + snapshots[n]+"/"+str(z)+"_densite_smooth2.fits"

                hdul = fits.open(data)
                field = hdul[0].data
                hdul.close()

                sum_ = np.sum(field[X0-25:X0+25,Y0-25:Y0+25,Z0-1:Z0+1],axis=2)


                if n > 0:
                    im = axes.imshow(sum_-sum_lcdm, origin="lower",vmin=-0.5,vmax = 0.5,cmap="bwr")
                else : 
                    im = axes.imshow(sum_, origin="lower",cmap="viridis")

                ims.append(im)

                if not k-1 > 5 : axes.xaxis.set_visible(False)  
                if not (k-1) % 3 == 0 : axes.yaxis.set_visible(False)  

                #axes.set_xlim(0,50)
                #axes.set_ylim(0,50)

                axes.set_xlabel(r"$\rm X [Mpc / h]$")
                axes.set_ylabel(r"$\rm Y [Mpc / h]$")


            im_lcdm = axs[0].images[0]  
            im_diff = axs[1].images[0]     

            divider_left = make_axes_locatable(axs[0])
            cax_left = divider_left.append_axes("left", size="5%", pad=0.1)
            plt.colorbar(im_lcdm, cax=cax_left)
            cax_left.yaxis.set_ticks_position('left')
            cax_left.yaxis.set_label_position('left')

            divider_right = make_axes_locatable(axs[-1])
            cax_right = divider_right.append_axes("right", size="5%", pad=0.1)
            plt.colorbar(im_diff, cax=cax_right)


            plt.tight_layout()
            plt.savefig(f"field_diff_{["min","max"][minmax]}_{Redshifts[i]}_{Xmax}.pdf")
            plt.clf()