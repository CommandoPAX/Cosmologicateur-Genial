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
    def __init__ (self, path, name = "lcdm", index = 2,tout = True):
        self.path = path
        self.json_path = ""

        self.index = index
        self.info = "info_"+"0"*(5-len(str(self.index//10)))+str(self.index)
        self.output = "output_"+"0"*(5-len(str(self.index//10)))+str(self.index)

        self.name = name

        self.args = {}        

        # Obliger de devoir charger l'intégralité des données pour calculer sigma_8 malheureusement
        
        self.data = yt.load(self.path+"/"+self.output+"/"+self.info+".txt")
        self.ad = self.data.all_data()

        self.df = self.ad.to_dataframe(["particle_position_x","particle_position_y","particle_position_z","particle_mass"])

        for i in self.df.columns :
            if not i in self.args : self.args[i] = self.df[i].to_numpy()
                
        #################################################################################
        # Constante qui sont utilisé pour Sigma_8, a terme devront nécessiter une récup depuis la base de données
        
        for root, dir, files in os.walk(self.path) : 
            for filename in files : 
                json_name= rf"{self.index}_PAR_.*\.json"
                if re.search(json_name, filename) : 
                    with open(f"{self.path}/{filename}", 'r') as f : 
                        self.Parameters = json.load(f) 
        
        self.BoxSize  = 500
        print(self.Parameters["namelist"])
        self.grid     = 2**self.Parameters["namelist"]["amr_params"]["levelmin"]                   #grid size
        self.MAS      = 'CIC'                   #Cloud-in-Cell
        self.verbose  = False   #whether print information on the progress   
        self.omega_m = self.Parameters["omega_m"]
        self.N=2**self.Parameters["namelist"]["amr_params"]["levelmin"]                 #grid size                           
        self.kmin=2*np.pi/self.BoxSize
        
                
        ###################################################################################

        self.Power_Spectrum()
        try :
             self.Halo()
        except :
             pass
         
        self.CST = {}
        self.CST["name"] = self.name

        self.Calc_Sigma_8()
        self.Create_Json()

    def __getitem__ (self, x):
        return self.args[x]

    def Create_Json(self) : 
        self.json_path = f"./RESULT/{self.name}_constants.json"
        with open(self.json_path, "w") as outf : 
            json.dump(self.CST, outf, indent=4, separators=(", ", ": "), sort_keys=True, skipkeys=True, ensure_ascii=False) 

    def Power_Spectrum(self) :
        for root,dirs,files in os.walk(self.path): 
                for file in files:
                    if "POW" in file and ".txt" in file and str(self.index)+"_" in file :
                        self.path_pow = file
                        break

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
    
    #############################################################################
    # Functions needed to calculate sigma_8
    
    def ifft(self, field):
        field[0,self.N//2,self.N//2]=0
        return np.fft.irfftn(np.fft.ifftshift(field.transpose()[:-1,:-1],axes=(0,1)),(self.N,self.N,self.N) )
    
    def fft(self, f_field): # fast fourier transform
        field=np.zeros((self.N//2+1,self.N+1,self.N+1),dtype=complex)
        field[:,:-1,:-1]=np.fft.fftshift(np.fft.rfftn(f_field),axes=(0,1)).transpose()
        field[:,-1],field[:,:,-1]=field[:,0],field[:,:,0]
        return field
    
    def W(self, grid):
        l = self.BoxSize/self.N
        k1,k2,k3=grid[0][self.N//2:]*l/2/np.pi,grid[1]*l/2/np.pi,grid[2]*l/2/np.pi
        return (np.sinc(k1)*np.sinc(k2)*np.sinc(k3))**2
    
    def WW(self, grid):
        R=8
        Rkmod=np.sqrt(grid[0][self.N//2:]**2+grid[1]**2+grid[2]**2)*R
        return 3 * (np.sin(Rkmod) - Rkmod*np.cos(Rkmod))/Rkmod**3

    def sigma8(self, field,k_grid):
        fdelta = self.fft(field)

        Wcic = self.W(k_grid)
        Wth = self.WW(k_grid)

        fdelta = fdelta * Wth / Wcic

        delta = self.ifft(fdelta)
        return np.sqrt(np.mean(delta**2))
        
    def Calc_Sigma_8(self) :
        # This function only works if yt has loaded all data
        
        k=np.linspace(-(self.N//2)*self.kmin,self.N//2*self.kmin,self.N+1,dtype=np.float64)
        k_grid=np.array(np.meshgrid(k,k,k,sparse=True,indexing='ij'),dtype=object) 
        
        delta = np.zeros((self.N,self.N,self.N), dtype=np.float32)
        
        pos = self.ad['particle_position'].astype(np.float32)*self.BoxSize
        MASL.MA(pos.astype(np.float32), delta, self.BoxSize, self.MAS, verbose=self.verbose)
        delta /= np.mean(delta, dtype=np.float64);  delta -= 1.0
        
        self.CST["sigma_8"] = self.sigma8(delta,k_grid)
        
        S_8 = self.CST["sigma_8"] * np.sqrt(self.omega_m / 0.3)
        self.CST["S_8"] = S_8

        print(S_8,self.CST["sigma_8"])

def Particle_mass (Simu):
    
    PPM = yt.ParticlePlot(Simu.data, 'particle_position_x', 'particle_position_y','particle_mass')
    PPM.set_unit('particle_mass', 'Msun')
    PPM.set_zlim(1e12,1e15)
    PPM.annotate_scale()

def Velocity (Simu) :
    
    VEL = yt.ParticlePlot(Simu.data, 'particle_velocity_x', 'particle_velocity_y','particle_mass')
    VEL.set_unit('particle_velocity_x', 'km/s')
    VEL.set_unit('particle_velocity_y', 'km/s')
    VEL.set_unit('particle_mass', 'Msun')

def Potential (Simu):
    
    POT = yt.SlicePlot(Simu.data, "z",('gravity', 'Potential'),center=[0.5, 0.5, 0.3])

def Plot_sigma_8 ():
    plt.clf()
    plt.figure()
    i = 0
    for root, dirs, files in os.walk("./RESULT/"):
        for file in files :
            if "_constants.json" in file :
                with open("./RESULT/"+file, "r") as f :
                    para = json.load(f)
                    S8 = para["S_8"]
                    plt.scatter([S8],[i], label = para["name"])
                    i = i+1
    axes = plt.gca()
    axes.set_xlim(0.2,0.6)
    axes.get_yaxis().set_visible(False)
    plt.xlabel(r'$S_8 \equiv \sigma_8 \sqrt{\Omega_m / 0.3}$', fontsize = 14)
    plt.legend()

    plt.savefig("./RESULT/S_8.png")

def Plot_Pow (Simu1, Ref, labelname : str = "Ratio") :
    plt.loglog(Ref["k"],Simu1,label=labelname) #plot measure from N-body
    axes = plt.gca()

    plt.xlabel("k [h/Mpc]")
    plt.ylabel(r"P(k) [$(Mpc/h)^3$]")
    plt.legend()
                    

if __name__ == "__main__" :
    
    cosmology.setCosmology('planck18')
    
    Path_lcdm = "./RESULT/2024-03-12 20:07:08 - LCDM" 

    lcdm = Simulation(Path_lcdm,name="lcdm",index = 2,tout = False)
    #noms = ["WDM500","WDM4000"]
    Pow = [0,0,0]
    Diviser_Pow(lcdm,lcdm)
    try : 
        for root, dirs, files in os.walk("./RESULT/"):

            for dir in dirs :

                nom = dir.split(" ")[-1]

                if "PGN" in nom and not "WDM" in nom : #in noms or noms == "":
                    simu2 = Simulation("./RESULT/"+dir,name=nom,index = 3,tout = False)
                    
                    """PowerSpectrum(lcdm)
                    PowerSpectrum(simu2)

                    plt.savefig("./RESULT/lcdm + "+nom+".png")
                    plt.clf()
                    if nom == "WDM3PGN1000" :
                        Pow[0] = simu2["Pk0"]/lcdm["Pk0"]
                    if nom == "PGN1000" :
                        Pow[1] = simu2["Pk0"]/lcdm["Pk0"]
                    if nom == "WDM3" :
                        Pow[2] = simu2["Pk0"]/lcdm["Pk0"]"""
                    
                    # Diviser_Pow(simu2,lcdm,ls="--")

        """plt.title ("z = 1")
        plt.savefig("./RESULT/WDM-nul-1.png")
                    #plt.clf()"""
    except : 
        pass 

Plot_sigma_8()

"""plt.clf()  
    print(Pow[0])
    print('#######################################################################')
    print(Pow[1])
    print('#######################################################################')
    print(Pow[2])
    print(Pow[0]-Pow[1]-Pow[2])
    Plot_Pow(Pow[0], lcdm, labelname = "WDM3PGN1000/lcdm")
    Plot_Pow(Pow[1] + Pow[2], lcdm, labelname = "(WDM3+PGN1000)/lcdm")
    Plot_Pow(-Pow[0] + (Pow[1] + Pow[2]), lcdm, labelname= "Différence entre les ")
    plt.savefig("./RESULT/RATIOOOOOO")
    Plot_sigma_8()   """             
