#    _      _                   _                 _                     _____ ______ _   _ _____          _      
#   | |    ( )                 | |               | |                   / ____|  ____| \ | |_   _|   /\   | |     
#   | |    |/  __ _ _ __   __ _| |_   _ ___  __ _| |_ ___ _   _ _ __  | |  __| |__  |  \| | | |    /  \  | |     
#   | |       / _` | '_ \ / _` | | | | / __|/ _` | __/ _ \ | | | '__| | | |_ |  __| | . ` | | |   / /\ \ | |     
#   | |____  | (_| | | | | (_| | | |_| \__ \ (_| | ||  __/ |_| | |    | |__| | |____| |\  |_| |_ / ____ \| |____ 
#   |______|  \__,_|_| |_|\__,_|_|\__, |___/\__,_|\__\___|\__,_|_|     \_____|______|_| \_|_____/_/    \_\______|
#                                  __/ |                                                                         
#                                 |___/      



import matplotlib.pyplot as plt
import yt
import numpy as np
from HaloStats import halo_MF
from colossus.cosmology import cosmology
from colossus.lss import mass_function
import os
import json
import matplotlib
import re 
import MAS_library as MASL

#matplotlib.rcParams.update({'font.size': 20})
#matplotlib.rcParams['text.usetex'] = True
#matplotlib.rcParams["figure.facecolor"]='w'
#matplotlib.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'

plt.tight_layout()

#aout = 0.09, 0.2, 0.333, 0.5, 1

def z(a):
    return 1/a - 1

index_redhift = {
    "2" : 10,
    "3" : 4,
    "4" : 2,
    "5" : 1,
    "6" : 0
}

class Simulation ():
    def __init__ (self, path, name = "LCDM", index = 2,tout = True):
        self.path = path
        self.json_path = ""

        self.index = index
        self.info = "info_"+"0"*(5-len(str(self.index//10)))+str(self.index)
        self.output = "output_"+"0"*(5-len(str(self.index//10)))+str(self.index)

        self.name = name

        self.args = {}        
        
        self.CST = {}
        self.CST["name"] = self.name


        self.Power_Spectrum()


    def __getitem__ (self, x):
        return self.args[x]


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

    def Convertir_index(self) : 

        # Convertit self.index en la vraie valeure
        # i = 2 -> avant dernier
        # i = 3 -> dernier

        para_file  = self.path + "/1_PAR_" + self.name +".json"
        with open(para_file, "r") as f :
            para = json.load(f)
            
        n_output = len(para["namelist"]["output_params"]["aout"])

        if self.index == 2 :
            self.index = para["namelist"]["output_params"]["aout"][n_output-2]
        elif self.index == 3 :
            self.index = para["namelist"]["output_params"]["aout"][n_output-1]

def PowerSpectrum (Simu, Class = False) :

    plt.loglog(Simu["k"],Simu["Pk0"],label=Simu.name) #plot measure from N-body
    plt.legend()

    if Class :
        toL=np.transpose(np.loadtxt("CLASS.dat"))
        plt.loglog(toL[0],toL[1],linestyle="dotted",label='CLASS') #plot lienar CLASS
        toL=np.transpose(np.loadtxt("CLASS_NL.dat"))
        plt.loglog(toL[0],toL[1],linestyle="dashdot",label='CLASS_NL') #plot non-linear CLASS from HaloFit
    
    plt.xlabel("k [h/Mpc]")
    axes = plt.gca()
    axes.set_xlim(2e-2,0.9)
    plt.ylabel(r"P(k) [$(Mpc/h)^3$]")
    plt.legend()


def Diviser_Pow (Simu1, Simu2,ls="") :
    if not ls :plt.loglog(Simu1["k"],Simu1["Pk0"]/Simu2["Pk0"],label="Ratio "+Simu1.name+" / "+Simu2.name) #plot measure from N-body
    else : plt.loglog(Simu1["k"],Simu1["Pk0"]/Simu2["Pk0"],label="Ratio "+Simu1.name+" / "+Simu2.name,ls=ls)
    axes = plt.gca()
    axes.set_xlim(2e-2,0.9)
    axes.set_ylim(0.6,1.2)
    plt.xlabel("k [h/Mpc]")
    plt.ylabel(r"P(k) [$(Mpc/h)^3$]")
    plt.legend()


def Plot_Pow (Simu1, Ref, labelname : str = "Ratio", linetype : str ="solid",color="", axes = plt) :
    if color == "" :axes.loglog(Ref["k"][8:160],Simu1,label=labelname, ls = linetype) #plot measure from N-body
    else :axes.loglog(Ref["k"][8:160],Simu1,label=labelname, ls = linetype,color=color) #plot measure from N-body


def trouver_simus (name, exclu = "", eq = True):

    result = []

    try : 
        for root, dirs, files in os.walk("./RESULT/"):

            for dir in dirs :

                nom = dir.split(" ")[-1]
                #print(nom, name, nom in name, nom in name and (not exclu in name or exclu == ""))

                if  name == nom and (not exclu in nom or exclu == "") and eq:
                    result.append("./RESULT/"+dir)

                if  name in nom and (not exclu in nom or exclu == "") and not eq:
                    result.append("./RESULT/"+dir)

    except : 
        pass           

    return result       

def superposer (fnl, wdm, path_lcdm = "./RESULT/2024-03-12 20:07:08 - LCDM"):

    lcdm = Simulation(path_lcdm,name="LCDM",index = 3,tout = False)
    lcdm2 = Simulation(path_lcdm,name="LCDM",index = 2,tout = False)

    simus = trouver_simus("WDM"+str(wdm)+"PGN"+str(fnl))
    print(simus)

    for s in simus :
        if "WDM"+str(wdm)+"P" in s.split(" ")[-1] : 
            pw2 = Simulation(s,name=s.split(" ")[-1], tout=False,index=2)
            pw3 = Simulation(s,name=s.split(" ")[-1], tout=False,index=3)

    simus = trouver_simus("WDM"+str(wdm), exclu="PGN")
    print(simus)
    for s in simus :
        if "WDM"+str(wdm) == s.split(" ")[-1] :
            w2 = Simulation(s,name=s.split(" ")[-1], tout=False,index=2)
            w3 = Simulation(s,name=s.split(" ")[-1], tout=False,index=3)

    simus = trouver_simus("PGN"+str(fnl), exclu="WDM")
    print(simus)
    for s in simus :
        p2 = Simulation(s,name=s.split(" ")[-1], tout=False,index=2)
        p3 = Simulation(s,name=s.split(" ")[-1], tout=False,index=3)

    pow_p2 = (p2["Pk0"]/lcdm2["Pk0"])[8:160]
    pow_w2 = (w2["Pk0"]/lcdm2["Pk0"])[8:160]
    pow_pw2 = (pw2["Pk0"]/lcdm2["Pk0"])[8:160]

    pow_p3 = (p3["Pk0"]/lcdm["Pk0"])[8:160]
    pow_w3 = (w3["Pk0"]/lcdm["Pk0"])[8:160]
    pow_pw3 = (pw3["Pk0"]/lcdm["Pk0"])[8:160]

    plt.clf()


    fig, axes = plt.subplots(nrows=2,figsize=(8,8))
    axes = axes.flatten()
    axes[0].set_xlim(0.1,2)
    #axes[0].set_ylim(pow_pw2[160],1.05)

    axes[1].set_xlim(0.1,2)

    axes[1].set_xlabel("k [h/Mpc]")
    axes[0].set_ylabel(r"P(k) [$(Mpc/h)^3$]")
    axes[1].set_ylabel(r"P(k) [$(Mpc/h)^3$]")


    div = pow_pw2 / (pow_p2 * pow_w2)

    #axes[1].set_ylim(np.min(div[8:160]),np.max(div[8:160]))

    Plot_Pow(pow_pw2, lcdm2, labelname = "(m = "+str(wdm)+" ev, fnl = "+str(fnl)+")/lcdm",axes=axes[0])
    Plot_Pow(pow_p2 * pow_w2, lcdm2, labelname = "((m = "+str(wdm)+" fnl = 0) * (m = 0 fnl = "+str(fnl)+")/lcdm",axes=axes[0])
    Plot_Pow(div, lcdm, labelname= "Ratio", axes = axes[1],color="black")
    Plot_Pow(pow_p2, lcdm2, labelname="(fnl = "+str(fnl)+")/lcdm", linetype='dotted',axes=axes[0],color="red")
    Plot_Pow(pow_w2, lcdm2, labelname="(m = "+str(wdm)+")/lcdm", linetype='dotted',color="black",axes=axes[0])

    axes[0].legend()
    axes[1].legend()

    plt.title("z = 1")

    plt.savefig("./RESULT/Superposition wdm"+str(wdm)+"fnl"+str(fnl)+" - z=1 .png")

    plt.clf()


    fig, axes = plt.subplots(nrows=2,figsize=(8,8))
    axes = axes.flatten()
    axes[0].set_xlim(0.1,2)
    #axes[0].set_ylim(pow_pw3[160],1.05)
    axes[1].set_xlim(0.1,2)

    axes[1].set_xlabel("k [h/Mpc]")
    axes[0].set_ylabel(r"P(k) [$(Mpc/h)^3$]")
    axes[1].set_ylabel(r"P(k) [$(Mpc/h)^3$]")



    plt.title("z = 0")

    div = pow_pw3 / (pow_p3 * pow_w3)

    #axes[1].set_ylim(np.min(div[8:160]),np.max(div[8:160]))


    Plot_Pow(pow_pw3, lcdm, labelname = "(m = "+str(wdm)+" ev, fnl = "+str(fnl)+"/lcdm",axes=axes[0])
    Plot_Pow((pow_p3 * pow_w3), lcdm, labelname = "((m = "+str(wdm)+" fnl = 0) * (m = 0 fnl = "+str(fnl)+")/lcdm",axes=axes[0])
    Plot_Pow(div, lcdm, labelname= "Ratio", axes = axes[1],color="black")
    Plot_Pow(pow_p3, lcdm, labelname="(fnl = "+str(fnl)+")/lcdm", linetype='dotted',axes=axes[0],color="red")
    Plot_Pow(pow_w3, lcdm, labelname="(m = "+str(wdm)+")/lcdm", linetype='dotted',color="black",axes=axes[0])
    axes[0].legend()
    axes[1].legend()
    plt.savefig("./RESULT/Superposition wdm"+str(wdm)+"fnl"+str(fnl)+" - z=0 .png")


def ractions (wdm = 100000000000000, fnl = True) :

    print("\n\nWDM "+str(wdm)+"\n\n")

    plt.title("z = 0")

    for root, dirs, files in os.walk("./RESULT"):
        for dir in dirs :
            if "WDM"+str(wdm) in dir and not "WDM"+str(wdm)+"0" in dir and "r" in dir :
                print(dir)
                s = Simulation("./RESULT/"+dir, name=dir.split(" ")[-1],index = 3, tout=False)
                Diviser_Pow(s,lcdm)

    plt.savefig("./RESULT/z = 0 - ractions - wdm"+str(wdm)+".png")

    plt.clf()

    plt.title("z = 1")

    for root, dirs, files in os.walk("./RESULT"):
        for dir in dirs :
            if "WDM"+str(wdm) in dir  and not "WDM"+str(wdm)+"0" in dir and "r" in dir :
                print(dir)
                s = Simulation("./RESULT/"+dir, name=dir.split(" ")[-1],index = 2, tout=False)
                Diviser_Pow(s,lcdm2)

    plt.savefig("./RESULT/z = 1 - ractions - wdm"+str(wdm)+".png")
    plt.clf()


if __name__ == "__main__" :
    

    cosmology.setCosmology('planck18')
    
    Path_lcdm = "./RESULT/2024-03-12 20:07:08 - LCDM" 
    #superposer(fnl=1000,wdm=300)
    #superposer(fnl=-1000,wdm=300)
    lcdm = Simulation(Path_lcdm, name="LCDM",index = 3,tout=False)
    lcdm2 = Simulation(Path_lcdm, name="LCDM",index = 2,tout=False)


    """for i in range(1,6):
        superposer(fnl=1000,wdm=i*100,path_lcdm=Path_lcdm)
        superposer(fnl=-1000,wdm=i*100,path_lcdm=Path_lcdm)"""
    

    ractions(10, True)    
    ractions(100, True)    
    ractions(200, True)    

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