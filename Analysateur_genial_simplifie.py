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
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams["figure.facecolor"]='w'
matplotlib.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'

plt.tight_layout()

#matplotlib.rc('font', family="serif",size=24)

#aout = 0.09, 0.2, 0.333, 0.5, 1

def z(a):
    return 1/a - 1

def ltok(x):
    return 2*np.pi/x

index_redhift = {
    "0" : 32,
    "1" : 3,
    "2" : 1,
    "3" : 0.25,
    "4" : 0
}

class Simulation ():
    def __init__ (self, name = "LCDM"):
        self.json_path = ""

        self.name = name

        print(name)
        
        self.CST = {}
        self.CST["name"] = self.name

        self.POW = []

        self.k0 = []
        self.k1 = []
        self.k2 = []
        self.k3 = []
        self.k4 = []

        self.P0 = []
        self.P1 = []
        self.P2 = []
        self.P3 = []       
        self.P4 = []

        self.k=[self.k0,self.k1,self.k2,self.k3,self.k4]
        self.P=[self.P0,self.P1,self.P2,self.P3,self.P4]

        for i in range(0,15):
            try :
                print(i)
                fichier = open("./POW/"+str(i)+"_POW_"+name+".txt","r")
                self.POW.append("./POW/"+str(i)+"_POW_"+name+".txt")
                fichier.close()
            except:
                print("err")
                break


        self.Power_Spectrum()


    def Power_Spectrum(self) :
        
        for i in range(5) :
            fichier = open(self.POW[i],"r")

            while 1:
                l = fichier.readline()
                if l == "" : break
                l = l.split(" ")

                self.k[i].append(float(l[0]))
                self.P[i].append(float(l[1].replace("\n","")))
                
            fichier.close()
 
        for i in range(len(self.P)) :
            self.P[i] = np.array(self.P[i])
            self.k[i] = np.array(self.k[i])
       

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
    
    lcdm = Simulation(name="benchM")


    WDM500 = Simulation(name="G_ViVi")
    WDM500fm500 = Simulation(name="NG_ViVi")
    WDM500f500 = Simulation(name="NG_Fminus500_ViVi")
    fm500 = Simulation(name="NG_F500")
    f500 = Simulation(name="NG_Fminus500")

    PNG_BM = Simulation(name="NG_F500")
    PNG_1 = Simulation(name="NG_F500_1")
    PNG_2 = Simulation(name="NG_F500_2")
    PNG_3 = Simulation(name="NG_F500_3")
    PNG_4 = Simulation(name="NG_F500_4")


    #Plot_sigma_8()
    

    n = 0
    
    Redshifts = {
        0 : 32,
        1 : 3,
        2 : 1,
        3 : 0.25,
        4 : 0
    }


    for i in [1]:


        #plt.subplot(2,2,n+1)

        axes = plt.gca()

        axes.title.set_text (f"z = {Redshifts[i]}")

        plt.loglog(lcdm.k[i], lcdm.P[i]/lcdm.P[i], label=r"$\Lambda{\rm CDM}$")
        #plt.loglog(WDM500.k[i], WDM500.P[i]/lcdm.P[i], label =r"$m_{\rm WDM} = 500 eV$", ls="--", color="green")
        #plt.loglog(WDM500f500.k[i], WDM500f500.P[i]/lcdm.P[i], label =r"$fnl = 500 \& WDM$", color="fuchsia",ls="--")
        #plt.loglog(WDM500fm500.k[i], WDM500fm500.P[i]/lcdm.P[i], label =r"$fnl = -500 \& WDM$", color="orange",ls="--")
        #plt.loglog(f500.k[i], f500.P[i]/lcdm.P[i], label =r"fnl = 500",color="fuchsia")
        #plt.loglog(fm500.k[i], fm500.P[i]/lcdm.P[i], label =r"$fnl = -500$", color="orange")
        
        plt.loglog(f500.k[i], f500.P[i]/lcdm.P[i], label =r"fnl = 500 1")
        plt.loglog(PNG_1.k[5], PNG_1.P[5]/lcdm.P[i], label =r"fnl = 500 2")
        plt.loglog(PNG_2.k[5], PNG_2.P[5]/lcdm.P[i], label =r"fnl = 500 3")
        plt.loglog(PNG_3.k[5], PNG_3.P[5]/lcdm.P[i], label =r"fnl = 500 4")
        plt.loglog(PNG_4.k[5], PNG_4.P[5]/lcdm.P[i], label =r"fnl = 500 5")

        
        n +=1

        if i == 0 or i == 2 : axes.set_ylabel(r"P/P$_{\rm CDM}$")

        if i == 0 : plt.legend()


    #axes.set_xlim(1e-2,3)


    
    print("test")
    plt.savefig("/home/fcastillo/Cosmologicateur-Genial/PS.pdf")
    

    #superposer(fnl=-1000,wdm=100)
    #superposer(fnl=1000,wdm=100)
    #superposer(fnl=-1000,wdm=200)

 