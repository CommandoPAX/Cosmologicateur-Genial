##### Le configurateur génial
##### Automatise la configuration de Ramses et monofonic


import os
import sys
import getopt

def generer_monofonic (argv) :

    opts, args = getopt.getopt(argv,"n:l:m:f:k:s:",["ngrid =","Lbox =","wdmass =","fnl =","kmin =","sigma ="])

    ramses = open("./ramses/namelist/ramses.nml","r")


    """
    n -> ngrid
    l -> taille de la boite (Mpc / h)
    m -> masse de la matière noire chaude (= 0 => LCDM)
    f -> fnl
    k -> kmin
    s -> sigma
    (f, k, s = 0 => pas de non gaussianités)
    """


    ngrid = 0
    taille = 0
    wdmass = 0  
    fnl = 0
    kmin = 0
    sigma = 0  


    for name, value in opts:
        if name in ["-n", "--ngrid"]:
            ngrid = value
        if name in ["-l", "--Lbox"]:
            taille = value
        if name in ["-w", "--wdmass"]:
            wdmass = value
        if name in ["-f", "--fnl"]:
            fnl = value
        if name in ["-k", "--kmin"]:
            kmin = value
        if name in ["-s", "--sigma"]:
            sigma = value

    SETUP = """
    [setup]

    """

    SETUP += "GridRes         = "+str(2**int(ngrid))+"""      # number of grid cells per linear dimension for calculations 
                            #   = particles for sc initial load
    """

    SETUP += "BoxLength       = "+str(taille)+"      # length of the box in Mpc/h"


    SETUP += """
    zstart          = 32.0     # starting redshift

    LPTorder        = 2        # order of the LPT to be used (1,2 or 3)

    DoBaryons       = no       # also do baryon ICs?
    DoBaryonVrel    = no       # if doing baryons, incl. also relative velocity to linear order?

    DoFixing        = no      # do mode fixing à la Angulo&Pontzen (https://arxiv.org/abs/1603.05253)
    DoInversion     = no       # invert phases (for paired simulations)

    ParticleLoad    = sc       # particle load, can be 'sc' (1x), 'bcc' (2x) or 'fcc' (4x) 
                            # (increases number of particles by given factor!), or 'glass'

    ## if `ParticleLoad = glass' then specify here where to load the glass distribution from
    # GlassFileName   = glass128.hdf5
    # GlassTiles      = 1

    #########################################################################################
    """


    COSMOLOGY = """
    [cosmology]
    ## transfer = ... specifies the Einstein-Boltzmann plugin module

    ParameterSet    = none # Planck2018EE+BAO+SN  # specify a pre-defined parameter set, or set to 'none' and set manually below
    """


    if wdmass !=0 :

        COSMOLOGY += """
    Omega_c         = 0.0
    N_ncdm          = 1     
    Omega_ncdm      = 0.26377    
    m_ncdm          = 1000      

    Omega_b         = 0.0494

    Omega_L         = 0.6842

    WDMmass="""+str(wdmass)+"""
    """    


    else: 
        COSMOLOGY += """
    Omega_c         = 0.26067  
    N_ncdm          = 0       
    Omega_ncdm      = 0.0     
    m_ncdm          = 0       
    Omega_b         = 0.04897 
    Omega_L         = 0.6889 

    """


    COSMOLOGY += """
    H0              = 67.66
    n_s             = 0.9665
    #sigma_8         = 0.8102
    A_s             = 2.1052e-9  # can use A_s instead of sigma_8 when using CLASS 
    norm            = 1.0
    Tcmb            = 2.7255
    k_p             = 0.05
    N_ur            = 2.046
    m_nu1           = 0.06
    m_nu2           = 0.0
    m_nu3           = 0.0
    w_0             = -1.0  # not supported yet!
    w_a             = 0.0   # not supported yet!
    transfer        = CLASS          
    ztarget         = 2.5             # target redshift for CLASS module, output at ztarget will be back-scaled to zstart
    """


    if int(fnl) !=0 and int(sigma) !=0 and int(kmin) !=0:

        COSMOLOGY += "fnl = "+str(fnl)+"\n"
        COSMOLOGY += "sigma = "+str(sigma)+"\n"
        COSMOLOGY += "kmin = "+str(kmin)+"\n"



    COSMOLOGY+="""
    ######################################################################################### 
    """

    RANDOM = """
    [random]
    ## generator = ... specifies the random field generator plugin module

    ##> NGenIC compatible random number generator module compatible with V. Springel's original code
    ## (https://www.h-its.org/2014/11/05/ngenic-code/) as well as the 2LPT code by Pueblas&Scoccmiarro
    ## (https://cosmo.nyu.edu/roman/2LPT/)
    generator      = NGENIC
    seed           = 12345

    ##> The PANPHASIA generator uses a plugin based on original code by A. Jenkins
    ## Warning: Before using this module, please make sure you read and agree to the distinct license
    ## requirements by registering on the website http://icc.dur.ac.uk/Panphasia.php

    # generator      = PANPHASIA
    # descriptor     = [Panph1,L10,(800,224,576),S9,CH1564365824,MXXL]
    # PanphasiaMinRootResolution = 512 # requires the white noise reallisation to be made at least at that resolution (default is 512)

    ##> The MUSIC1 multi-scale random number generator is provided for convenience
    ## warning: MUSIC1 generator is not MPI parallel (yet) (memory is needed for full field on each task)
    # generator      = MUSIC1  
    # seed[7]        = 12345
    # seed[8]        = 23456
    # seed[9]        = 34567

    # Add a possible constraint field here:
    # ConstraintFieldFile = initial_conditions.hdf5
    # ConstraintFieldName = ic_white_noise


    #########################################################################################"""
    EXECUTION = """
    [execution]
    # Specify the number of threads / task
    NumThreads      = 8
    """
    OUTPUT = """

    #########################################################################################
    [output]
    ## format = .... specifies the output plugin module

    ##> RAMSES / GRAFIC2 compatible format
    format	        = grafic2
    filename        = ics_ramses
    grafic_use_SPT  = no # if no then uses PPT, otherwise linear SPT

    ##> Gadget-2/3 'fortran unformatted binary'-style format
    # format          = gadget2
    # filename        = ics_gadget.dat
    # UseLongids      = false

    ##> Gadget-2/3 HDF5 format
    # format          = gadget_hdf5
    # filename        = ics_gadget.hdf5

    ##> Arepo HDF5 format (virtually identical to gadget_hdf5)
    # format          = AREPO
    # filename        = ics_arepo.hdf5

    ##> HACC compatible generic-io format
    # format          = genericio
    # filename        = ics_hacc

    ##> SWIFT compatible HDF5 format. Format broadly similar to gadget_hdf5 but in a single
    ##> file even when using MPI. No h-factors for position and masses and no sqrt(a)-factor for the velocities.
    ##> IDs are stored using 64-bits unless UseLongids is set to false.
    # format          = SWIFT
    # filename        = ics_swift.hdf5
    # UseLongids      = true

    ##> Generic HDF5 output format for testing or PT-based calculations
    # format          = generic
    # filename        = debug.hdf5
    # generic_out_eulerian = yes  # if yes then uses PPT for output
    """

    MONOFONIC = SETUP + COSMOLOGY + RANDOM + EXECUTION + OUTPUT


    monofonic = open("./monofonic_exp/config.conf","w")
    monofonic.write(MONOFONIC)
    monofonic.close()

    output_ramses = ""
    while 1 :
        ligne = ramses.readline()
        if ligne == "" : break
        for i in range(10): ligne = ligne.replace("  "," ")
        ligne = ligne.split("=")
        if ligne[0] == "levelmin" : ligne[1] = ngrid+"\n"
        if ligne[0] == "levelmax" : ligne[1] = str(int(ngrid)+6)+"\n"

        for k in ligne :
            output_ramses += k
            if not "\n" in k : 
                output_ramses+="="

        output_ramses = output_ramses.replace("initfile(1)='monofonic/","initfile(1)='monofonic_exp/")

    ramses.close()

    ramses = open("./ramses/namelist/ramses.nml","w")
    ramses.write(output_ramses)
    ramses.close()



def compiler ():

    os.system("rm -fr ./monofonic_exp/build")
    os.system("mkdir ./monofonic_exp/build")
    os.system("cd ./monofonic_exp/build ; cmake .." )
    os.system("cd ./monofonic_exp/build ; make")
    os.system("cd ./monofonic_exp; build/monofonIC config.conf")
    os.system("cp ./monofonic_exp/ics_ramses/ic_poscx ./monofonic_exp/ics_ramses/ic_deltab")
    os.system("cd ./ramses/bin; make clean")
    os.system("cd ./ramses/bin; make NDIM=3")
    os.system("cp ./ramses/bin/ramses3d ..")
    os.system("cp ./ramses/namelist/ramses.nml ..")

    #os.system("./ramses/bin/ramses3d ./ramses/namelist/ramses.nml")  le fichier bash le fait déjà

if __name__ == "__main__" :
    generer_monofonic(sys.argv[1:])
    compiler()