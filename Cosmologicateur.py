import matplotlib.pyplot as plt
import matplotlib

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
import os

cosmology.setCosmology('planck18')

ds=yt.load('../output_00002/info_00002.txt')

#plot
q = yt.ParticlePlot(ds, 'particle_position_x', 'particle_position_y','particle_mass')
q.set_unit('particle_mass', 'Msun')
#q.zoom(4)
#q.annotate_timestamp(corner='upper_left', time=True, redshift=False, draw_inset_box=True,time_format='t = {time:.1f}', time_unit='code_time')
q.annotate_scale()
q.show()
#q.save()

