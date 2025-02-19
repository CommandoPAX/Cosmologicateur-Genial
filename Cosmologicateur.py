#    _                                          _             _           _                     _____ ______ _   _ _____          _      
#   | |                                        | |           (_)         | |                   / ____|  ____| \ | |_   _|   /\   | |     
#   | |     ___    ___ ___  ___ _ __ ___   ___ | | ___   __ _ _  ___ __ _| |_ ___ _   _ _ __  | |  __| |__  |  \| | | |    /  \  | |     
#   | |    / _ \  / __/ _ \/ __| '_ ` _ \ / _ \| |/ _ \ / _` | |/ __/ _` | __/ _ \ | | | '__| | | |_ |  __| | . ` | | |   / /\ \ | |     
#   | |___|  __/ | (_| (_) \__ \ | | | | | (_) | | (_) | (_| | | (_| (_| | ||  __/ |_| | |    | |__| | |____| |\  |_| |_ / ____ \| |____ 
#   |______\___|  \___\___/|___/_| |_| |_|\___/|_|\___/ \__, |_|\___\__,_|\__\___|\__,_|_|     \_____|______|_| \_|_____/_/    \_\______|
#                                                        __/ |                                                                           
#

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
from Errorateur import LogError

# Utilisation (not currently functionnal):
# python Cosmologicateur.py -i entree -o sortie

Result_Path = "./RESULT" #Path where all results will be saved, default is Cosmologicateur-Genial/RESULT/
#yt.enable_parallelism()

def Copy_Mono_Config(path : str) : 

    os.system(f'cp ./monofonic_exp/config.conf "{path}/config.conf"')

def Predicted_Particle_Mass(DATA, index : int, path : str, SimuName : str) :

    if SimuName == "" : 
        output_ = f"{path}/{index}_PPM.png"
    else : 
        output_ = f"{path}/{index}_PPM_{SimuName}.png"
    #plot
    PPM = yt.ParticlePlot(DATA, 'particle_position_x', 'particle_position_y','particle_mass')
    PPM.set_unit('particle_mass', 'Msun')
    try :
        PPM.set_zlim(('particle_mass'),zmin=(5e12,"Msun"),zmax=(2e14,"Msun"))
        PPM.annotate_title(SimuName)
    except :
        pass
    #PPM.zoom(4)
    #PPM.annotate_timestamp(corner='upper_left', time=True, redshift=False, draw_inset_box=True,time_format='t = {time:.1f}', time_unit='code_time')
    PPM.annotate_scale()
    PPM.save(output_)


def Potential(DATA, index : int, path : str, SimuName : str) :
    
    if SimuName == "" : 
        output_ = f"{path}/{index}_POT.png"
    else : 
        output_ = f"{path}/{index}_POT_{SimuName}.png"
    POT = yt.SlicePlot(DATA, "z",('gravity', 'Potential'),center=[0.5, 0.5, 0.3])
    #POT.annotate_cell_edges()
    POT.save(output_)

def Velocity(DATA, index : int, path : str, SimuName : str) :
    
    if SimuName == "" : 
        output_ = f"{path}/{index}_VEL.png"
    else : 
        output_ = f"{path}/{index}_VEL_{SimuName}.png"
    VEL = yt.ParticlePlot(DATA, 'particle_velocity_x', 'particle_velocity_y','particle_mass')
    VEL.set_unit('particle_velocity_x', 'km/s')
    VEL.set_unit('particle_velocity_y', 'km/s')
    VEL.set_unit('particle_mass', 'Msun')
    VEL.save(output_)

def Power_Spectrum(DATA, index : int, path : str, SimuName : str) :

    # Define important parameters
    if SimuName == "" : 
        output_ = f"{path}/{index}_POW.png"
    else : 
        output_ = f"{path}/{index}_POW_{SimuName}.png"

    grid = 512    #grid size
    pBoxSize = DATA.domain_width.in_units('Mpc/h') #Mpc/h
    BoxSize = pBoxSize[0].value #Mpc/h
    print(BoxSize)
    Rayleigh_sampling = 1     #whether sampling the Rayleigh distribution for modes amplitudes
    threads = 1      #number of openmp threads
    verbose = True   #whether to print some information
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
        

    Pk = PKL.Pk(delta, BoxSize, axis, MAS, threads, verbose)
    k       = Pk.k3D
    Pk0     = Pk.Pk[:,0]

    try :
        output_spectrum = open(output_[:-4]+".txt","w")

        for i in range(len(k)):
            output_spectrum.write(str(k[i])+" "+str(Pk0[i])+"\n")

        output_spectrum.close()
    except :
        pass


    plt.clf()

    plt.loglog(k,Pk0,label="RAMSES") #plot measure from N-body

    toL=np.transpose(np.loadtxt("CLASS.dat"))
    plt.loglog(toL[0],toL[1],linestyle="dotted",label='CLASS') #plot lienar CLASS
    toL=np.transpose(np.loadtxt("CLASS_NL.dat"))
    plt.loglog(toL[0],toL[1],linestyle="dashdot",label='CLASS_NL') #plot non-linear CLASS from HaloFit
    plt.legend()
    plt.xlabel("k [h/Mpc]")
    plt.ylabel(r"P(k) [$(Mpc/h)^3$]")
    plt.savefig(output_)


def Halo(DATA, index, path : str, SimuName : str) : 

    try :
        if SimuName == "" : 
            output_ = f"{path}/{index}_HAL.png"
        else : 
            output_ = f"{path}/{index}_HAL_{SimuName}.png"
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

        os.system(f"cp ./halo_catalogs/info_0000{index}/info_0000{index}.0.h5 {path}")

        bin_centers, num_halos, err = halo_MF(haloM, log_M_min=log_M_min, log_M_max=log_M_max, delta_log_M=delta_log_M, boxsize = boxsize) #calculate halo mass function

        fig,ax = plt.subplots()
        ax.errorbar(bin_centers[num_halos!=0], num_halos[num_halos!=0], yerr=err[num_halos!=0], fmt='x', capthick=2, elinewidth=2,label='Measured') #plt HMF


        HMF_T=mass_function.massFunction(bin_centers, 0.0, mdef = 'fof', model = 'sheth99',q_out = 'dndlnM')      #get theoretical line for HMF,
        #check https://bdiemer.bitbucket.io/colossus/lss_mass_function.html for a list of model
        plt.loglog(bin_centers, HMF_T,label="Theo")

        ax.set_xscale('log')
        ax.set_yscale('log')

        plt.legend()
        plt.xlabel(r'Mass ($M_{\odot}$)]', fontsize = 14)
        plt.ylabel(r'$\frac{dN}{d\log M}$ ($Mpc^{-3}$)', fontsize = 14)
        plt.savefig(output_)
    except:
        pass


def Get_Simu_Info(DATA, index, path : str, SimuName : str) : #Not sure if there will be a different one for each dataset

    if SimuName == "" : 
        output_ = f"{path}/{index}_PAR.json"
    else : 
        output_ = f"{path}/{index}_PAR_{SimuName}.json"
    with open(output_, "w") as outf : 
        Simu_Info = DATA.parameters
        json.dump(Simu_Info, outf, indent=4, separators=(", ", ": "), sort_keys=True, skipkeys=True, ensure_ascii=False) 

def main(argv):
    global Result_Path
    
    POT = True
    VEL = True 
    SPE = True 
    HAL = False

    pre = "../../../data77/stahl/Scale/Nb/WDM/KF/"
    snapshots = ["benchM","NG_F500","G_m500","NG_F500_m500","NG_Fminus500","NG_Fminus500_m500"]

    n = int (sys.argv[1])
    i = int(sys.argv[2])

    name = snapshots[n]
    file_path = pre + name

    cosmology.setCosmology('planck18')
    
    os.system(f"mkdir {Result_Path}/{name}")

    Output_Path = f"{Result_Path}/{name}"
    

    print(f'---------------------------------{i}----------------------------------------')
    
    input_ = f"{file_path}/snapshot_00{i}.hdf5"

    units_override =  {"UnitLength_in_cm": 3.08568e21}

    ds=yt.load(input_,unit_base=units_override)
    
    #Predicted_Particle_Mass(ds, i, Output_Path, name)
    #Get_Simu_Info(ds, i, Output_Path, name)
    #if POT : Potential(ds, i, Output_Path, name)
    #if VEL : Velocity(ds, i, Output_Path, name)
    if SPE : Power_Spectrum(ds, i, Output_Path, name)
    if HAL : Halo(ds, i, Output_Path, name)
    #os.system(f'cp -r ./output_0000{i} "{Output_Path}/output_0000{i}"')

    
if __name__ == "__main__" :
    
    #pre = "../../../data77/stahl/Scale/Nb/WDM/KF/"
    #snapshots = ["NG_F500","G_m500","NG_F500_m500","NG_Fminus500_noScale","NG_Fminus500","NG_Fminus500_m500"]

    #pre = "../../../data77/stahl/Scale/Nb/WDM/ViVi/"
    #snapshots = ["G_ViVi","NG_ViVi","NG_Fminus500_ViVi"]):

    #snap = snapshots[0]
    #file_path = pre + snap

    main(sys.argv)
    
    #for snap in snapshots :
    #    
    #    file_path = pre + snap

    #    main(file_path, name = snap)
