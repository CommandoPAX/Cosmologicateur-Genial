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


pre = "/data77/stahl/Scale/Nb/WDM/KF/"
snapshots = ["benchM","NG_F500","G_m500","NG_F500_m500","NG_Fminus500","NG_Fminus500_m500", "G_ViVi", "NG_ViVi", "NG_Fminus500_ViVi"]
Result_Path = "/data100/fcastillo/RESULT/" 

for n in range(6):
    for i in range(1,5):
        print(n,i)
        file = h5py.File(pre+snapshots[n]+"/fof_subhalo_tab_00"+str(i)+".hdf5")
        print(file, file.keys())
        pos = file["Group"]["GroupPos"][:]
        file.close()
        output = open(Result_Path+snapshots[n]+"/"+str(i)+"_halos.txt","w")

        output.write("# X Y Z")
        pos = str(pos).replace("[","")
        pos = str(pos).replace("]","")
        output.write(pos)
        print(str(pos))

        output.close()