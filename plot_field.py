import numpy as np
import matplotlib.pyplot as plt
import matplotlib


matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams["figure.facecolor"]='w'
matplotlib.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'

snapshots = ["benchM","NG_F500","G_m500","NG_F500_m500","NG_Fminus500","NG_Fminus500_m500", "G_ViVi","NG_ViVi" , "NG_Fminus500_ViVi","NEDE", "NsPNG_EDE_F500","NsPNG_EDE_F1833"]

z= [15,12, 10, 8, 5,3,1,0.5,0.25,0]
indices_z = [5,6,8,9]

indices_hdm = [0,1,4,6,7,8,9,10,11]
#indices_hdm = [0,2,6]
#indices_hdm = [0,1,4,6,9]

for i in range(4):
    plt.figure()
    for n in indices_hdm: 
        if n <= 8 : redshifts = range(1,5)
        else : redshifts = indices_z
        for i in redshifts:
            z = redshifts[i]
            
