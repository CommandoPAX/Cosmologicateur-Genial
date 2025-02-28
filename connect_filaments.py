from skeletonnateur_hpc import*
import sys
import numpy as np

N = int(sys.argv[1])
i = int(sys.argv[2])

snapshots = ["benchM","NG_F500","G_m500","NG_F500_m500","NG_Fminus500","NG_Fminus500_m500"]
labels = ["LCDM", "fnl = -500", "m = 500 eV", "WDM & fnl = -500", "fnl = 500", "WDM & fnl = 500"]


squelette = Squelette_3d(f"/data100/fcastillo/RESULT/{snapshots[N]}/{i}_densite_0_c0.1.up.NDskl.S001.a.NDskl")

n = []
for p in squelette.Pointscrit :
    if int(p.type) == 3: n.append(p.nfil)

np.save(f"/data100/fcastillo/RESULT/{snapshots[N]}/{i}_densite_0_c0.1_connect_fil.txt", np.array(n))