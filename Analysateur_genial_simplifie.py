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
    def __init__ (self, name = "LCDM"):
        self.json_path = ""

        self.name = name

        print(name)
        
        self.CST = {}
        self.CST["name"] = self.name

        self.POW = []

        self.k2 = []
        self.k3 = []

        self.P2 = []
        self.P3 = []

        self.k=np.array([self.k2,self.k3])
        self.P=np.array([self.P2,self.P3])

        for i in range(1,6):
            try :
                fichier = open("./POW/"+str(i)+"_POW_"+name+".txt","r")
                self.POW.append("./POW/"+str(i)+"_POW_"+name+".txt")
                fichier.close()
            except:
                break


        self.Power_Spectrum()


    def Power_Spectrum(self) :
        
        for i in range(2) :
            fichier = open(self.POW[len(self.POW)-2+i])

            while 1:
                l = fichier.readline()
                if l == "" : break
                l = l.split(" ")

                if i == 0:
                    self.k2.append(float(l[0]))
                    self.P2.append(float(l[1].replace("\n","")))
                else :
                    self.k3.append(float(l[0]))
                    self.P3.append(float(l[1].replace("\n","")))

            fichier.close()

        self.P2 = np.array(self.P2)
        self.P3 = np.array(self.P3)


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

def superposer (fnl, wdm):

    global lcdm

    w = Simulation(name="WDM"+str(wdm))
    p = Simulation(name="PGN"+str(fnl))
    pw = Simulation(name="WDM"+str(wdm)+"PGN"+str(fnl))

    plt.loglog(lcdm.k2, w.P2/lcdm.P2, ls="dotted",color="red",label="WDM alone, m = "+str(wdm)+" ev")
    plt.loglog(lcdm.k2, p.P2/lcdm.P2, ls="dotted",color="black",label="NG alone, fnl = "+str(fnl))
    plt.loglog(lcdm.k2, p.P2*w.P2/(lcdm.P2**2), color="orange",label=r"Multiplication $P_w \times P_w$ "+str(wdm))
    plt.loglog(lcdm.k2, pw.P2/lcdm.P2, color="blue",label="WDM and NG, m = "+str(wdm)+" ev, fnl = "+str(fnl))

    axes = plt.gca()

    axes.set_xlabel("k [h/Mpc]")
    axes.set_ylabel(r"P(k) [$(Mpc/h)^3$]")

    plt.legend()

    plt.show()

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
    
    #superposer(fnl=1000,wdm=300)
    #superposer(fnl=-1000,wdm=300)
    lcdm = Simulation(name="LCDM")

    superposer(fnl=-1000,wdm=100)

    noms = []
    for root, dirs, files in os.walk("./POW"):
        for file in files :
            nom = file.replace("1_","")
            nom = nom.replace("2_","")
            nom = nom.replace("3_","")
            nom = nom.replace("4_","")
            nom = nom.replace("5_","")
            nom = nom.replace("6_","")

            nom = nom.replace("POW_","")
            nom = nom.replace(".txt","")

            if (not nom in noms) and ((not "r" in nom) and ("1000P" in nom) or "LCDM" in nom) : 
                noms.append(nom)

    for nom in noms :
        simu = Simulation(name=nom)
        plt.loglog(simu.k2, simu.P2/lcdm.P2,label=nom)


    plt.legend()
    plt.show()

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
