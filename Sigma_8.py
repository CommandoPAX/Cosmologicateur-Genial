# Handles data processing

import matplotlib
import matplotlib.pyplot as plt
import yt
import numpy as np
from colossus.cosmology import cosmology
import os
import json
import re 
import MAS_library as MASL

plt.rcParams['figure.figsize'] = [20, 7]


matplotlib.rcParams.update({'font.size': 20})

matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams["figure.facecolor"]='w'
matplotlib.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'

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
        
        """ Unused values
        self.A_s = 2.1064e-9
        self.n_s = 0.96822
        self.k_pivot = 0.05 # 1/Mpc
        self.omega_r= 0 #9.16714e-05
        self.omega_k=self.Parameters["omega_k"]
        self.z=50
        self.N_ncdm=0
        self.m_ncdm=0.0
        self.omega_ncdm= 0 #0.00140718
        self.omega_cdm=0.261008
        self.omega_b=self.Parameters["omega_b"]
        
        #self.h=0.67742 Don't know if that's H0 in PAR file but I will assume it is, divided by 100
        self.h = self.Parameters["H0"]/100
        
        self.axis     = 0                       #axis along which place RSD; not used here
        self.do_RSD   = False                   #dont do redshif-space distortions
        self.ptypes   = [1]                   #CDM + neutrinos
        """
                
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
        
if __name__ == "__main__" :
    
    cosmology.setCosmology('planck18')
    noms = ""#["WDM500","WDM4000"]
    if 1:#try : 
        for root, dirs, files in os.walk("./RESULT/"):

            for dir in dirs :

                nom = dir.split(" ")[-1]

                if True: #if nom in noms or noms == "":
                    plt.clf()
                    simu2 = Simulation("./RESULT/"+dir,name=nom,index = 3,tout = True)
    else:
        pass 