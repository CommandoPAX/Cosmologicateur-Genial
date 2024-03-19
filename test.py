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
DATA.all_data().to_dataframe(["particle_position_x","particle_position_y","particle_position_z","particle_mass"])

grid = 256    #grid size
pBoxSize = DATA.domain_width.in_units('Mpccm/h') #Mpc/h
BoxSize = pBoxSize[0].value #Mpc/h
Rayleigh_sampling = 1     #whether sampling the Rayleigh distribution for modes amplitudes
threads = 1      #number of openmp threads
verbose = False   #whether to print some information
axis = 0
MAS = 'CIC'

ad=DATA.all_data()
pos = ad['particle_position'].astype(np.float32)*BoxSize

# define 3D density fields
delta = np.zeros((grid,grid,grid), dtype=np.float32)

# construct 3D density field
MASL.MA(pos.astype(np.float32), delta, BoxSize, MAS, verbose=verbose)

# at this point, delta contains the effective number of particles in each voxel
# now compute overdensity and density constrast
delta /= np.mean(delta, dtype=np.float64);  delta -= 1.0

