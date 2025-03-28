import matplotlib.pyplot as plt
import yt
import numpy as np
from scipy import interpolate
import density_field_library as DFL
import Pk_library as PKL
import MAS_library as MASL
import mass_function_library as MFL
from yt_astro_analysis.halo_analysis import HaloCatalog
from HaloStats import halo_MF
from colossus.cosmology import cosmology
from colossus.lss import mass_function
import os, sys, getopt
import json, datetime 
import h5py
from astropy.io import fits


pre = "/data77/stahl/Scale/Nb/WDM/ViVi/"
snapshots = ["benchM","NG_F500","G_m500","NG_F500_m500","NG_Fminus500","NG_Fminus500_m500", "G_ViVi", "NG_ViVi", "NG_Fminus500_ViVi"]
Result_Path = "/data100/fcastillo/RESULT/" 

for n in range(6,9):
    for i in range(5):
        file = h5py.File(pre+snapshots[n]+"/fof_subhalo_tab_00"+str(i)+".hdf5")
        pos = file["Group"]["GroupPos"][:]
        file.close()
        header = fits.Header()
        header['COMMENT'] = 'Position halos'
        header['NAXIS'] = 3  
        header['NAXIS1'] = len(pos)  
        header['NAXIS2'] = len(pos) 
        header['NAXIS3'] = len(pos)

        hdu = fits.PrimaryHDU(data=pos, header=header)

        hdu.writeto(f"{Result_Path}/{i}_halos.fits", overwrite=True)
