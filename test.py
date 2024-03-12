import matplotlib.pyplot as plt
import yt
import numpy as np
from scipy import interpolate
import density_field_library as DFL
import Pk_library as PKL
import MAS_library as MASL
import mass_function_library as MFL
from yt.extensions.astro_analysis.halo_analysis import HaloCatalog
from HaloStats import halo_MF
from colossus.cosmology import cosmology
from colossus.lss import mass_function
import os, sys, getopt
import json, datetime 
from Errorateur import LogError

index = 2
DATA = yt.load("../test/output_00002/info_00002.txt")
pBoxSize = DATA.domain_width.in_units('Mpc/h') #Mpc/h
BoxSize = pBoxSize[0].value #Mpc/h
hc = HaloCatalog(data_ds=DATA, finder_method="hop") #Run halo Finder
hc.create()
ds = yt.load(f"./halo_catalogs/info_0000{index}/info_0000{index}.0.h5") #Get the file saved by hc.create
ad = ds.all_data()
# The halo mass
haloM=ad["halos", "particle_mass"]        
log_M_min=14.3 #minimal mass to plot HMF
log_M_max=15 #maximal mass to plot HMF
delta_log_M=0.1
boxsize=BoxSize/0.67 #factor h
try :
    os.system(f"cp ./halo_catalogs/info_0000{index}/info_0000{index}.0.h5 {output_[-4:]}.h0.5")
except:
    pass

bin_centers, num_halos, err = halo_MF(haloM, log_M_min=log_M_min, log_M_max=log_M_max, delta_log_M=delta_log_M, boxsize = boxsize) #calculate halo mass function
print(bin_centers, num_halos,err)
