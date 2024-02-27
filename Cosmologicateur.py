# Handles data processing

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

# Utilisation (not currently functionnal):
# python Cosmologicateur.py -i entree -o sortie

Result_Path = "./RESULT" #Path where all results will be saved, default is Cosmologicateur-Genial/RESULT/
Ramses_Path = "../ramses"
Data_Path = "../../Results/LCDM_256"

def Copy_Mono_Config(path : str) : 
    try : 
        os.system(f'cp ../monofonic_exp/config.conf "{path}/config.conf"')
    except Exception as e : 
        LogError("Copy_Mono_Config", e)
        print(e)

def Predicted_Particle_Mass(DATA, index : int, path : str) :
    try :
        output_ = f"{path}/{index}_PPM.png"
        #plot
        PPM = yt.ParticlePlot(DATA, 'particle_position_x', 'particle_position_y','particle_mass')
        PPM.set_unit('particle_mass', 'Msun')
        #PPM.zoom(4)
        #PPM.annotate_timestamp(corner='upper_left', time=True, redshift=False, draw_inset_box=True,time_format='t = {time:.1f}', time_unit='code_time')
        PPM.annotate_scale()
        PPM.save(output_)
    except Exception as e : 
        LogError("Predicted_Particle_Mass", e)
        print(e)

def Potential(DATA, index : int, path : str) : 
    try : 
        output_ = f"{path}/{index}_POT.png"
        POT = yt.SlicePlot(DATA, "z",('gravity', 'Potential'),center=[0.5, 0.5, 0.3])
        POT.annotate_cell_edges()
        POT.save(output_)
    except Exception as e : 
        LogError("Potential", e)
        print(e)

def Velocity(DATA, index : int, path : str) : 
    try : 
        output_ = f"{path}/{index}_VEL.png"
        VEL = yt.ParticlePlot(DATA, 'particle_velocity_x', 'particle_velocity_y','particle_mass')
        VEL.set_unit('particle_velocity_x', 'km/s')
        VEL.set_unit('particle_velocity_y', 'km/s')
        VEL.set_unit('particle_mass', 'Msun')
        VEL.save(output_)
    except Exception as e : 
        LogError("Velocity", e)
        print(e)

def Power_Spectrum(DATA, index : int, path : str) :
    try :  
        # Define important parameters
        output_ = f"{path}/{index}_POW.png"
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
        
        Pk = PKL.Pk(delta, BoxSize, axis, 'None', threads, verbose)
        k       = Pk.k3D
        Pk0     = Pk.Pk[:,0]

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
    except Exception as e : 
        LogError("Power_Spectrum", e)
        print(e)

def Halo(DATA, index, path : str) : 
    try :
        output_ = f"{path}/{index}_HAL.png"
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
    except Exception as e : 
        LogError("Halo", e)
        print(e)

def Get_Simu_Info(DATA, index, path : str) : #Not sure if there will be a different one for each dataset
    try : 
        output_ = f"{path}/{index}_PAR.json"
        with open(output_, "w") as outf : 
            Simu_Info = DATA.parameters
            json.dump(Simu_Info, outf, indent=4, separators=(", ", ": "), sort_keys=True, skipkeys=True, ensure_ascii=False)
    except Exception as e : 
        LogError("Get_Simu_Info", e)
        print(e) 

def main(argv):
    global Result_Path
    global Ramses_Path
    
    POT = False 
    VEL = False 
    SPE = False 
    HAL = False 
    
    try:
        # Predicted particle mass will be enabled by default for the pretty pictures
        # p for potential 
        # v for velocity 
        # s for power spectrum 
        # m for halo mass (no need for help function)
        opts, args = getopt.getopt(argv,"pvsma")
    except getopt.GetoptError:
        print ('test.py -p -v -s -m')
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-p"):
            POT = True
        elif opt in ("-v"):
            VEL = True
        elif opt in ("-s"):
            SPE = True 
        elif opt in ("-m"): 
            HAL = True
        elif opt in ("-a"): 
            POT = True
            VEL = True 
            SPE = True 
            HAL = True 
    
    cosmology.setCosmology('planck18')
    
    Output_Path = f"{Result_Path}/{str(datetime.datetime.now())[:-7]}"
    
    for i in range(1, 10)  : 
        plt.clf()
        print(f'---------------------------------{i}----------------------------------------')
        if 1:#try : #Will load files until they don't exist anymore
            input_ = f"../output_0000{i}/info_0000{i}.txt"
            ramses_input_ = input_
            ds=yt.load(input_)
            rds = yt.load(ramses_input_)
        else:#except : 
            print("File not found, breaking Thomas Delzant's legs")
            break 
        Predicted_Particle_Mass(ds, i, Output_Path)
        Get_Simu_Info(rds, i, Output_Path)
        if POT : Potential(ds, i, Output_Path)
        if VEL : Velocity(ds, i, Output_Path)
        if SPE : Power_Spectrum(rds, i, Output_Path)
        if HAL : Halo(rds, i, Output_Path)
        os.system(f'cp -r ../output_0000{i} "{Output_Path}/ouput_0000{i}"')
        os.system(f'cp -r {Ramses_Path}/output_0000{i} "{Output_Path}/ramses_output_0000{i}"')
    if POT and VEL and HAL and SPE : Copy_Mono_Config(Output_Path) 

def test():
    global Data_Path

    cosmology.setCosmology('planck18')

    for i in range(2):
        input_ = f"{Data_Path}/info_0000{i+1}.txt"
        ramses_input_ = f"{Data_Path}/rinfo_0000{i+1}.txt"
        ds=yt.load(input_)
        rds = yt.load(ramses_input_)



if __name__ == "__main__" :
    main(sys.argv[1:])
