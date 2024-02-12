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
import sys, getopt

# Utilisation :
# python Cosmologicateur.py -i entree -o sortie


def main(argv):

    inputfile = '../output_00002/info_00002.txt'
    outputfile = "Resultat.png"

    try:
      opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
    except getopt.GetoptError:
        print ('test.py -i <inputfile> -o <outputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('test.py -i <inputfile> -o <outputfile>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg


    monofonic = open("../monofonic/monofonic.conf","r")
    lignes = []
    gridres = 0
    sizebox = 0

    while 1 :
        ligne = monofonic.readline()
        for i in range(10):
            ligne = ligne.replace ("  ", " ")
        if ligne =="" : break
        ligne = ligne.split(" ")
        if ligne[0] == "GridRes" : gridres = ligne[2]
        if ligne[0] == "BoxLength" : sizebox = ligne[2]

    outputfile = outputfile[:-4]+"-"+str(gridres)+outputfile[-4:]
    outputfile = outputfile[:-4]+"-"+str(sizebox)+" Mpc"+outputfile[-4:]

    cosmology.setCosmology('planck18')

    ds=yt.load(inputfile)

    #plot
    q = yt.ParticlePlot(ds, 'particle_position_x', 'particle_position_y','particle_mass')
    q.set_unit('particle_mass', 'Msun')
    #q.zoom(4)
    #q.annotate_timestamp(corner='upper_left', time=True, redshift=False, draw_inset_box=True,time_format='t = {time:.1f}', time_unit='code_time')
    q.annotate_scale()
    #q.show()
    q.save(outputfile)

if __name__ == "__main__" :
    main(sys.argv[1:])