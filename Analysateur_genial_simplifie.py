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

matplotlib.rc('font', family="serif",size=24)

#aout = 0.09, 0.2, 0.333, 0.5, 1

def z(a):
    return 1/a - 1

def ltok(x):
    return 2*np.pi/x

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

def Plot_sigma_8 (index = 3, name="WDM",exclu="PGN"):
    plt.clf()
    i = 0
    liste_fichiers = ["lcdm",'WDM100',"WDM200","WDM300","WDM400","WDM500","PGN500","PGN1000","PGN-500","PGN-1000"]
    colors = ["#000000","#00FF00","#00DD00","#00BB00","#009900","#007700","#CC0000","#FF0000","#0000CC","#0000FF"]
    for root, dirs, files in os.walk("../json/"):
        for j in range(len(liste_fichiers)):
            for file in files :
                if file == "2_"+liste_fichiers[j]+"_constants.json" :
                    try : 
                        with open("../json/"+file, "r") as f :
                            para = json.load(f)
                            S8 = para["S_8"]
                            if True : #if (name in para["name"] and not exclu in para["name"]) or "cdm" in para["name"] :

                                if "PGN" in para ["name"] :
                                    legende = para["name"].replace("PGN",r"$f_{\rm nl}^0 = $")
                                    print(legende)
                                elif "WDM" in para ["name"] :
                                    legende = para["name"].replace("WDM",r"WDM, $m = $")
                                    legende+=r" ev"
                                elif "lcdm" in para["name"]:
                                    legende = r"$\Lambda$CDM"
                                else :
                                    print(para["name"])
                                    legende = para["name"]
                                plt.scatter([S8],[i], label = legende,s=150,color=colors[j])
                                i = i+1
                                if "cdm" in para["name"]:
                                    plt.axvline(x = S8,ls="--",color="black")
                    except:
                        pass

    axes = plt.gca()
    axes.set_xlim(0.6,1)
    axes.get_yaxis().set_visible(False)
    plt.xlabel(r'$S_8 \equiv \sigma_8 \sqrt{\Omega_m / 0.3}$', fontsize = 32)
    handles, labels = axes.get_legend_handles_labels()
    axes.legend(handles[::-1], labels[::-1], loc='upper left')
    #plt.legend()
    plt.show()


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

    #plt.title("z = 1")

    # {\rm WDM} m=100 {\rm ev}, {\rm NG}, {\rm NG} \times {\rm WDM} 

    plt.loglog(lcdm.k2, w.P2/lcdm.P2, ls="dotted",color="red",label=r"${\rm WDM}$, m=100 ${\rm ev}$")
    plt.loglog(lcdm.k2, p.P2/lcdm.P2, ls="dotted",color="black",label=r"${\rm NG}$, fnl = "+str(fnl))
    plt.loglog(lcdm.k2, p.P2*w.P2/(lcdm.P2**2), color="orange",label=r"${\rm NG} \times {\rm WDM}$ ")
    plt.loglog(lcdm.k2, pw.P2/lcdm.P2, color="blue",label=r"${\rm NG + WDM}$")

    axes = plt.gca()

    secax = axes.secondary_xaxis('top', functions=(ltok, ltok))
    secax.set_xlabel(r'$L{\rm [Mpc}/h]$') 

    plt.axvline(x = 2*np.pi/(500/0.67)*256, color = 'k')

    axes.set_xlabel("k [h/Mpc]")
    axes.set_ylabel(r"P / P$_{\Lambda {\rm CDM}}$")

    plt.legend()

    plt.show()

    """plt.title("z = 0")

    plt.loglog(lcdm.k3, w.P3/lcdm.P3, ls="dotted",color="red",label=r"${\rm WDM}$, m=100 ${\rm ev}$")
    plt.loglog(lcdm.k3, p.P3/lcdm.P3, ls="dotted",color="black",label=r"${\rm NG}$, fnl = "+str(fnl))
    plt.loglog(lcdm.k3, p.P3*w.P3/(lcdm.P3**2), color="orange",label=r"${\rm NG} \times {\rm WDM}$")
    plt.loglog(lcdm.k3, pw.P3/lcdm.P3, color="blue",label=r"${\rm NG + WDM}$")

    axes = plt.gca()

    plt.axvline(x = 2*np.pi/(500/0.67)*512, color = 'k')

    secax = axes.secondary_xaxis('top', functions=(ltok, ltok))
    secax.set_xlabel(r'$L{\rm [Mpc}/h]$') 

    axes.set_xlabel("k [h/Mpc]")
    axes.set_ylabel(r"$P(k) / P_{lcdm} (k)$")

    plt.legend()

    plt.show()"""


if __name__ == "__main__" :
    

    cosmology.setCosmology('planck18')
    
    lcdm = Simulation(name="LCDM")

    Plot_sigma_8()

    """WDM100 = Simulation(name="WDM100")
    WDM500 = Simulation(name="WDM500")
    WDM1000 = Simulation(name="WDM1000")
    WDM200 = Simulation(name="WDM200")
    WDM300 = Simulation(name="WDM300")

    Plot_sigma_8()
    
    plt.loglog(lcdm.k2, lcdm.P2/lcdm.P2, label=r"$\Lambda{\rm CDM}$")
    plt.loglog(WDM100.k2, WDM100.P2/lcdm.P2, label =r"$m = 100$")
    plt.loglog(WDM200.k2, WDM200.P2/lcdm.P2, label =r"$m = 200$")
    plt.loglog(WDM400.k2, WDM400.P2/lcdm.P2, label =r"$m = 300$")
    plt.loglog(WDM500.k2, WDM500.P2/lcdm.P2, label =r"$m = 500$")
    plt.loglog(WDM1000.k2, WDM1000.P2/lcdm.P2, label =r"$m = 1000$")
    

    superposer(fnl=-1000,wdm=100)
    superposer(fnl=1000,wdm=100)
    superposer(fnl=-1000,wdm=200)

   
plt.clf()  
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
