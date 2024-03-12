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

ds = yt.load("../test/output_00002/info_00002.txt")
print(ds.all_data())
plot = yt.ParticlePlot(ds, 'particle_position_x', 'particle_position_y','particle_mass')
plot.save("test.pdf")
os.system("evince test.pdf")

fichier = open("../fields","r")
while 1:
    l = fichier.readline()
    if "’, ‘" in l : print(l)
    if l == "" : break
