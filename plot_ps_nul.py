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
matplotlib.rcParams.update({'font.size': 20})

toL=np.transpose(np.loadtxt("CLASS_NL.dat"))
plt.loglog(toL[0],toL[1],label=r"$\Lambda{\rm CDM}$") #plot non-linear CLASS from HaloFit

axes = plt.gca()

axes.set_xlabel(r"${\rm Wavenumber}~k~[h/{\rm Mpc]}$")
axes.set_ylabel(r"$P(k) ~[ ({\rm Mpc} / h)^3]$")

plt.legend()
plt.tight_layout()
plt.savefig("/home/fcastillo/Cosmologicateur-Genial/pow_lcdm.pdf")