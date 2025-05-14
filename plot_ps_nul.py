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


fichier = open("/home/fcastillo/Cosmologicateur-Genial/4_POW_benchM.txt","r")

k = []
P = []

while 1:
    l = fichier.readline()
    if l == "" : break
    l = l.split(" ")

    k.append(float(l[0]))
    P.append(float(l[1].replace("\n","")))
    
fichier.close()

plt.plot(k,P)
axes = plt.gca()

axes.set_xlabel(r"${\rm Wavenumber}~k{\rm~[h/Mpc]}$")
axes.set_ylabel(r"$P(k) {\rm (Mpc / h)}^3$")

plt.savefig("/home/fcastillo/Cosmologicateur-Genial/pow_lcdm.pdf")