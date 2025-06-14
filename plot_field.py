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

snapshots = ["benchM","NG_F500","G_m500","NG_F500_m500","NG_Fminus500","NG_Fminus500_m500", "G_ViVi","NG_ViVi" , "NG_Fminus500_ViVi","NEDE", "NsPNG_EDE_F500","NsPNG_EDE_F1833","NBM"]
labels = [r"$\Lambda$CDM", r"$f_{\rm NL} = -500$", "m = 500 eV", "WDM & fnl = -500", r"$f_{\rm NL} = 500$", "WDM & fnl = 500", r"$m_{\rm WDM} = 10~{\rm eV}, f_{\rm WDM} = 2~\%$", r"$f_{\rm NL} = -500~\&~{\rm mixed~DM}$", r"$f_{\rm NL} = 500~\&~{\rm mixed~DM}$", r"${\rm EDE}$", r"$f_{\rm NL} = -300~\&~{\rm EDE}$", r"$f_{\rm NL} = -1100~\&~{\rm EDE}$",r"$\Lambda{\rm CDM}$"]

z= [15,12, 10, 8, 5,3,1,0.5,0.25,0]
indices_z = [5,6,8,9]

indices_hdm = [0,1,4,6,7,8]
#indices_hdm = [0,2,6]
#indices_hdm = [0,1,4,6,9]

Redshifts = [3,1,0.25,0]
pre = "/data100/fcastillo/RESULT/"



for BM in range(2):
    if BM == 0 : 
        indices_hdm = [0,1,4,6,7,8]
        nrow = 2
        ncol = 3
    else :
        indices_hdm = [12,9,10,11]
        nrow = 2
        ncol = 2
    for minmax in range(2) :
        ims = []

        for Xmax in [50]:
            for i in range(4):

                outer = gridspec.GridSpec(nrows=nrow, ncols=ncol)

                plt.figure(figsize=(9,6))


                axs = []
                for row in range(nrow):
                    for col in range(ncol):
                        inner = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=1, subplot_spec=outer[row, col], hspace=0)
                        axs += [plt.subplot(cell) for cell in inner]


                k = 0

                if BM == 0 : data = pre + snapshots[0]+"/"+str(i+1)+"_densite_smooth2.fits"
                else : data = pre + snapshots[12]+"/"+str(redshifts[i])+"_densite_smooth2.fits"
        
                hdul = fits.open(data)
                lcdm = hdul[0].data
                hdul.close()

                if minmax == 0:
                    max_index_flat = np.argmax(lcdm)
                    max_position = np.unravel_index(max_index_flat, lcdm.shape)
                else:
                    P = 40
                    extrema = np.load(f"/data100/fcastillo/RESULT/extrema/extrema_{0}_{i}_{2}.txt.npy")
                    minimums = extrema[extrema[:, 3] == 3][:, 0:3] % 512

                    # coords filtrées loin du bord
                    margin_xy = 25
                    margin_z = 1

                    # Indices valides (dans le cube 512x512x512)
                    valid = (
                        (minimums[:, 0] >= margin_xy) & (minimums[:, 0] <= 512 - margin_xy) &
                        (minimums[:, 1] >= margin_xy) & (minimums[:, 1] <= 512 - margin_xy) &
                        (minimums[:, 2] >= margin_z)  & (minimums[:, 2] <= 512 - margin_z)
                    )
                    minimums = minimums[valid]

                    Xk, Yk, Zk = minimums[:, 0].astype(int), minimums[:, 1].astype(int), minimums[:, 2].astype(int)

                    vals = lcdm[Xk, Yk, Zk]
                    seuil_haut_k = np.percentile(vals, 5)
                    seuil_bas_k = np.percentile(vals, 1)

                    mask = (vals >= seuil_bas_k) & (vals <= seuil_haut_k)
                    points_k = minimums[mask]

                    max_position = tuple(points_k[0].astype(int))

                X0, Y0, Z0 = max_position
                print(X0, Y0, Z0)


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


                    data = pre + snapshots[n]+"/"+str(z)+"_densite_smooth2.fits"

                    hdul = fits.open(data)
                    field = hdul[0].data
                    hdul.close()

                    sum_ = np.sum(field[X0-25:X0+25,Y0-25:Y0+25,Z0-1:Z0+1],axis=2)


                    if n > 0 and n!=12:
                        axes.title.set_text(labels[n]+r"$ - \Lambda{\rm CDM}$")
                        im = axes.imshow(sum_-sum_lcdm, origin="lower",vmin=-0.5,vmax = 0.5,cmap="bwr")
                    else : 
                        axes.title.set_text(r"$\Lambda{\rm CDM}$")
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
                plt.colorbar(im_diff, cax=cax_right,pad = 0.05)

                plt.subplots_adjust(left=0.15)
                plt.tight_layout()
                plt.savefig(f"field_diff_{["max","min"][minmax]}_{Redshifts[i]}_{Xmax}_{["","EDE"][BM]}.pdf")
                plt.clf()