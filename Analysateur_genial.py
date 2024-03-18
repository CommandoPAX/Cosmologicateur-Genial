# Handles data processing

import matplotlib.pyplot as plt
import yt
import numpy as np
import mass_function_library as MFL
from HaloStats import halo_MF
from colossus.cosmology import cosmology
from colossus.lss import mass_function
import os
import json
import re 

class Simulation ():
    def __init__ (self, path, name = "lcdm", index = 2):
        self.path = path
        self.json_path = ""

        self.index = index
        self.info = "info_"+"0"*(5-len(str(self.index//10)))+str(self.index)
        self.output = "output_"+"0"*(5-len(str(self.index//10)))+str(self.index)

        self.name = name

        self.data = yt.load(self.path+"/"+self.output+"/"+self.info+".txt")
        self.ad = self.data.all_data()

        self.df = self.ad.to_dataframe(["particle_position_x","particle_position_y","particle_position_z","particle_mass"])

        self.args = {}        

        for i in self.df.columns :
            if not i in self.args : self.args[i] = self.df[i].to_numpy()


        self.Power_Spectrum()
        try :
             self.Halo()
        except :
             pass
         
        self.CST = {}

        self.Calc_Sigma_8()
        self.Create_Json()

    def __getitem__ (self, x):
        return self.args[x]

    def Power_Spectrum(self) :

        path = ""
        for root,dirs,files in os.walk(self.path): 
                for file in files:
                    if "POW" in file and ".txt" in file and str(self.index)+"_" in file :
                        self.path_pow = file
                        break
        print(self.path_pow)
        fichier = open(self.path+"/"+self.path_pow,"r")
        k = []
        Pk0 = []

        while 1:
            l = fichier.readline()
            if l == "" : break
            l = l.split(" ")
            k.append(float(l[0]))
            Pk0.append(float(l[1].replace("\n","")))

        self.k = np.array(k)
        self.Pk0 = np.array(Pk0)

        self.args["k"] = self.k
        self.args["Pk0"] = self.Pk0

    def Halo(self) : 

        for root,dirs,files in os.walk(self.path): 
                for file in files:

                    if ".0.h5" in file :
                        path_halo = file
                        break        

        ds = yt.load(path_halo) #Get the file saved by hc.create
        ad = ds.all_data()
        # The halo mass
        haloM=ad["halos", "particle_mass"]      

        self.halo = haloM.to_numpy()
        self.args["halos"] = self.halo
        
    def Create_Json(self) : 
        self.json_path = f"{self.path}/constants.json"
        with open(self.json_path, "w") as outf : 
            json.dump(self.CST, outf, indent=4, separators=(", ", ": "), sort_keys=True, skipkeys=True, ensure_ascii=False) 
            
    def Calc_Sigma_8(self) :
        # compute the value of sigma_8
        sigma_8 = MFL.sigma(self.args["k"], self.args["Pk0"], 8.0)
        self.CST["sigma_8"] = sigma_8
        for root, dir, files in os.walk(self.path) : 
            for filename in files : 
                if re.search(r"3_PAR_.*\.json", filename) : 
                    with open(f"{self.path}/{filename}", 'r') as f : 
                        Parameters = json.load(f) 
        Omega_m = Parameters["omega_m"]
        S_8 = sigma_8 * np.sqrt(Omega_m / 0.3)
        self.CST["S_8"] = S_8
        self.args["S_8"] = S_8
        self.args["sigma_8"] = sigma_8
                    
        
def PowerSpectrum (Simu, Class = False) :

    plt.loglog(Simu["k"],Simu["Pk0"],label=Simu.name) #plot measure from N-body
    plt.legend()

    if Class :
        toL=np.transpose(np.loadtxt("CLASS.dat"))
        plt.loglog(toL[0],toL[1],linestyle="dotted",label='CLASS') #plot lienar CLASS
        toL=np.transpose(np.loadtxt("CLASS_NL.dat"))
        plt.loglog(toL[0],toL[1],linestyle="dashdot",label='CLASS_NL') #plot non-linear CLASS from HaloFit
        plt.legend()
        plt.xlabel("k [h/Mpc]")
        plt.ylabel(r"P(k) [$(Mpc/h)^3$]")

def Diviser_Pow (Simu1, Simu2) :
    plt.loglog(Simu1["k"],Simu1["Pk0"]/Simu2["Pk0"],label="Ratio "+Simu1.name+" / "+Simu2.name) #plot measure from N-body
    plt.legend()
            
def Halo (Simu, log_M_min = 14.3, log_M_max=15,delta_log_M=0.1) :
     

    pBoxSize = Simu.data.domain_width.in_units('Mpc/h') #Mpc/h
    BoxSize = pBoxSize[0].value #Mpc/h

    log_M_min=14.3 #minimal mass to plot HMF
    log_M_max=15 #maximal mass to plot HMF
    delta_log_M=0.1
    boxsize=BoxSize/0.67 #factor h

    bin_centers, num_halos, err = halo_MF(Simu["halos"], log_M_min=log_M_min, log_M_max=log_M_max, delta_log_M=delta_log_M, boxsize = boxsize) #calculate halo mass function

    fig,ax = plt.subplots()
    ax.errorbar(bin_centers[num_halos!=0], num_halos[num_halos!=0], yerr=err[num_halos!=0], fmt='x', capthick=2, elinewidth=2,label='Measured_'+str(Simu.name)) #plt HMF


    HMF_T=mass_function.massFunction(bin_centers, 0.0, mdef = 'fof', model = 'sheth99',q_out = 'dndlnM')      #get theoretical line for HMF,
    plt.loglog(bin_centers, HMF_T,label="Theo_"+str(Simu.name))

    ax.set_xscale('log')
    ax.set_yscale('log')

    plt.legend()
    plt.xlabel(r'Mass ($M_{\odot}$)]', fontsize = 14)
    plt.ylabel(r'$\frac{dN}{d\log M}$ ($Mpc^{-3}$)', fontsize = 14)

def Particle_mass (Simu):
    
    PPM = yt.ParticlePlot(Simu.data, 'particle_position_x', 'particle_position_y','particle_mass')
    PPM.set_unit('particle_mass', 'Msun')
    PPM.annotate_scale()

def Velocity (Simu) :
    
    VEL = yt.ParticlePlot(Simu.data, 'particle_velocity_x', 'particle_velocity_y','particle_mass')
    VEL.set_unit('particle_velocity_x', 'km/s')
    VEL.set_unit('particle_velocity_y', 'km/s')
    VEL.set_unit('particle_mass', 'Msun')

def Potential (Simu):
    
    POT = yt.SlicePlot(Simu.data, "z",('gravity', 'Potential'),center=[0.5, 0.5, 0.3])

def Ouvrir_sigma_8 ():
    plt.figure()
    i = 0
    for root, dirs, files in os.walk("./RESULT/"):
        for file in files :
            if re.search(r"constants\.json", file) :
                with open("./RESULT/"+file, "r") as f :
                    para = json.load(f)
                    plt.plot((i, para["S_8"]), label = para["name"])
                    plt.legend()
                    
                    

if __name__ == "__main__" :

    #sigma_8 ()
    
    cosmology.setCosmology('planck18')

    Path_lcdm = "./RESULT/2024-03-12 20:07:08 - LCDM" 

    lcdm = Simulation(Path_lcdm,name="lcdm",index = 3)
    for root, dirs, files in os.walk("./RESULT/"):

        for dir in dirs :

            nom = dir.split(" ")[-1]

            if not "LCDM" in nom :
                simu2 = Simulation("./RESULT/"+dir,name=nom,index = 3)
                
                PowerSpectrum(lcdm)
                PowerSpectrum(simu2)

                plt.savefig("lcdm + "+nom+".png")
                plt.clf()

                Diviser_Pow(lcdm, simu2)
                plt.savefig("lcdm + "+nom+"-d.png")
                plt.clf()
