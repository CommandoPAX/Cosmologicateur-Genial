from skeletonnateur_hpc import*
import sys
import numpy as np

n = int(sys.argv[1])
i = int(sys.argv[2])

snapshots = ["benchM","NG_F500","G_m500","NG_F500_m500","NG_Fminus500","NG_Fminus500_m500"]
labels = ["LCDM", "fnl = -500", "m = 500 eV", "WDM & fnl = -500", "fnl = 500", "WDM & fnl = 500"]


squelette = Squelette_2d(f"/data100/fcastillo/RESULT/{snapshots[n]}/{i}_densite_0_c0.1.up.NDskl.S001.a.NDskl")

longueurs = []
for fil in squelette.liste_filaments :
    fil.len()
    longueurs.append(fil.l)

hist = np.histogram(np.array(longueurs), density= True, range = [0, 10], bins=20)
np.save(f"/data100/fcastillo/RESULT/{snapshots[n]}/{i}_densite_0_c0.1_len_fil.txt", hist[0])