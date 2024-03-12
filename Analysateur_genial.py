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

class Simulation ():
    def __init__ (self, path, output_path = "", name = "lcdm",index = 2):
        self.path = path
        self.output_path = output_path

        self.index = index
        self.info = "info_"+"0"*(4-len(str(self.index//10)))+str(self.index)
        self.output = "output_"+"0"*(4-len(str(self.index//10)))+str(self.index)

        self.name = name

        self.data = yt.load(self.path)

        if self.output_path =="":
            self.output_Path = f"{self.path}/{str(datetime.datetime.now())[:-7]} - {self.name}
    
    def Copy_Mono_Config(self) : 

        os.system(f'cp ./monofonic_exp/config.conf {self.output_path}/config.conf')

    def Predicted_Particle_Mass(self, save = True) :

        self.PPM = yt.ParticlePlot(self.data, 'particle_position_x', 'particle_position_y','particle_mass')
        self.PPM.set_unit('particle_mass', 'Msun')
        self.PPM.annotate_scale()
        
        if save:
            output_ = f"{self.output_path}/PPM_{self.name}.png"
            self.PPM.save(output_)

        return self.PPM


    def Potential(self, save = True) :
        
        self.POT = yt.SlicePlot(self.data, "z",('gravity', 'Potential'),center=[0.5, 0.5, 0.3])
        self.POT.annotate_cell_edges()

        if save :
            output_ = self.output_path+"/POT_"+self.name+".png"
            self.POT.save(output_)

        return self.POT

    def Velocity(self, save = True) :
        
        self.VEL = yt.ParticlePlot(self.DATA, 'particle_velocity_x', 'particle_velocity_y','particle_mass')
        self.VEL.set_unit('particle_velocity_x', 'km/s')
        self.VEL.set_unit('particle_velocity_y', 'km/s')
        self.VEL.set_unit('particle_mass', 'Msun')
        
        if save :
            output_ = self.output_path+"/VEL_"+self.name+".png"
            self.VEL.save(output_)

        return self.VEL

    def Power_Spectrum(self, save = False, Class = False) :

        # Define important parameters

        grid = 256    #grid size
        pBoxSize = self.data.domain_width.in_units('Mpccm/h') #Mpc/h
        BoxSize = pBoxSize[0].value #Mpc/h
        Rayleigh_sampling = 1     #whether sampling the Rayleigh distribution for modes amplitudes
        threads = 1      #number of openmp threads
        verbose = False   #whether to print some information
        axis = 0
        MAS = 'CIC'
        
        ad=self.data.all_data()
        pos = ad['particle_position'].astype(np.float32)*BoxSize

        # define 3D density fields
        delta = np.zeros((grid,grid,grid), dtype=np.float32)

        # construct 3D density field
        MASL.MA(pos.astype(np.float32), delta, BoxSize, MAS, verbose=verbose)

        # at this point, delta contains the effective number of particles in each voxel
        # now compute overdensity and density constrast
        delta /= np.mean(delta, dtype=np.float64)
        delta -= 1.0
        
        Pk = PKL.Pk(delta, BoxSize, axis, 'None', threads, verbose)
        k       = Pk.k3D
        Pk0     = Pk.Pk[:,0]

        plt.clf()

        plt.loglog(k,Pk0,label="RAMSES") #plot measure from N-body

        if Class :
            toL=np.transpose(np.loadtxt("CLASS.dat"))
            plt.loglog(toL[0],toL[1],linestyle="dotted",label='CLASS') #plot lienar CLASS
            toL=np.transpose(np.loadtxt("CLASS_NL.dat"))
            plt.loglog(toL[0],toL[1],linestyle="dashdot",label='CLASS_NL') #plot non-linear CLASS from HaloFit
            
        plt.legend()
        plt.xlabel("k [h/Mpc]")
        plt.ylabel(r"P(k) [$(Mpc/h)^3$]")

        if save :
            output_ = self.output_path+"/POW_"+self.name+".png"
            plt.savefig(output_)


    def Halo(self, save = True) : 

        try :
            pBoxSize = self.data.domain_width.in_units('Mpc/h') #Mpc/h
            BoxSize = pBoxSize[0].value #Mpc/h
            os.system ("mkdir "+self.path+"/halo_catalogs")
            hc = HaloCatalog(data_ds=self.data, finder_method="hop",output_dir = self.path+"./halo_catalogs") #Run halo Finder
            hc.create()

            ds = yt.load(self.path+"/halo_catalogs/"+self.info+"/"+self.info+".0.h5") #Get the file saved by hc.create
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

            if save :

                output_ = self.output_path+"/HAL_"+self.name+".png"
                plt.savefig(output_)
        except:
            pass


    def Get_Simu_Info(self, save = False) : #Not sure if there will be a different one for each dataset

        output_ = self.output_path+"/PAR_"+self.name+".png"

        with open(output_, "w") as outf : 
            Simu_Info = self.data.parameters
            json.dump(Simu_Info, outf, indent=4, separators=(", ", ": "), sort_keys=True, skipkeys=True, ensure_ascii=False) 

def main(argv):
    global Result_Path
    
    POT = False 
    VEL = False 
    SPE = False 
    HAL = False 
    name = ""

    # Predicted particle mass will be enabled by default for the pretty pictures
    # p for potential 
    # v for velocity 
    # s for power spectrum 
    # m for halo mass (no need for help function)
    opts, args = getopt.getopt(argv,"pvsman:")

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
        elif opt in ("-n") :
            name = str(arg)
    
    cosmology.setCosmology('planck18')

    Result_Path = "./Cosmologicateur-Genial/RESULT/output_00001/info_00001.txt"

    if name =="" : 
        Output_Path = f"{Result_Path}/{str(datetime.datetime.now())[:-7]}"
    else : 
        Output_Path = f"{Result_Path}/{str(datetime.datetime.now())[:-7]} - {name}"
    
    lcdm = Simulation(Result_Path,name="lcdm",output_path=Output_Path)

    lcdm.Predicted_Particle_Mass()
    lcdm.Get_Simu_Info(False)
    if POT : lcdm.Potential(False)
    if VEL : lcdm.Velocity(False)
    if SPE : lcdm.Power_Spectrum(False)
    if HAL : lcdm.Halo(False)

    #os.system(f'cp -r {Ramses_Path}/output_0000{i} "{Output_Path}/ramses_output_0000{i}"')
    #lcdm.Copy_Mono_Config(Output_Path) 

if __name__ == "__main__" :
    main(sys.argv[1:])
