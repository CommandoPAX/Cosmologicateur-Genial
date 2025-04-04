from skeletonnateur_hpc import*
import sys
import numpy as np

n = int(sys.argv[1])
i = int(sys.argv[2])

#snapshots = ["benchM","NG_F500","G_m500","NG_F500_m500","NG_Fminus500","NG_Fminus500_m500", "G_ViVi", "NG_ViVi", "NG_Fminus500_ViVi"]
#labels = ["LCDM", "fnl = -500", "m = 500 eV", "WDM & fnl = -500", "fnl = 500", "WDM & fnl = 500"]

snapshots = ["NEDE","NsPNG_F500","NsPNG_F1000","NsPNG_F1833","NsPNG_EDE_F500","NsPNG_EDE_F1000","NsPNG_EDE_F1833"]

z= [15,12, 10, 8, 5,3,1,0.5,0.25,0]
indices_z = [5,6,8,9]



squelette = Squelette_3d(f"/data100/fcastillo/RESULT/{snapshots[n]}/{i}_densite_smooth2_c0.1.up.NDskl.S001.a.NDskl")

longueurs = []
for fil in squelette.liste_filaments :
    fil.len()
    longueurs.append(fil.l)

np.save(f"/data100/fcastillo/RESULT/{snapshots[n]}/{i}_densite_smooth2_c0.1_len_fil.txt", np.array(longueurs))
