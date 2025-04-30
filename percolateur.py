from skeletonnateur_hpc import*
import sys
import numpy as np
from astropy.io import fits

from scipy.spatial import KDTree

def identify_clusters(points, linking_length):
    tree = KDTree(points)
    N = len(points)
    visited = np.zeros(N, dtype=bool)
    clusters = []

    for i in range(N):
        if not visited[i]:
            cluster = []
            stack = [i]
            visited[i] = True

            while stack:
                idx = stack.pop()
                cluster.append(idx)
                neighbors = tree.query_ball_point(points[idx], r=linking_length)

                for neighbor in neighbors:
                    if not visited[neighbor]:
                        visited[neighbor] = True
                        stack.append(neighbor)

            clusters.append(cluster)

    return clusters

n = int(sys.argv[1])
i = int(sys.argv[2])

snapshots = ["benchM","NG_F500","G_m500","NG_F500_m500","NG_Fminus500","NG_Fminus500_m500", "G_ViVi","NG_ViVi" , "NG_Fminus500_ViVi"]

z= [15,12, 10, 8, 5,3,1,0.5,0.25,0]
indices_z = [5,6,8,9]

snapshots = ["benchM","NG_F500","G_m500","NG_F500_m500","NG_Fminus500","NG_Fminus500_m500", "G_ViVi","NG_ViVi" , "NG_Fminus500_ViVi","NEDE", "NsPNG_EDE_F500","NsPNG_EDE_F1833"]


R = 2
P = 5

result = np.load(f"/data100/fcastillo/RESULT/extrema/extrema_{n}_{i}_{R}.txt.npy")

pre = "/data100/fcastillo/RESULT/"

data = pre + snapshots[n]+"/"+str(i)+"_densite_smooth2.fits"

hdul = fits.open(data)
field = hdul[0].data
hdul.close()



X, Y, Z = result[:, 0].astype(int), result[:, 1].astype(int), result[:, 2].astype(int)


seuil_haut = np.percentile(field[X,Y,Z], P)

indices = field[X, Y, Z] >= seuil_haut


result = result[indices]
xyz = result[result[:,3]==0][:,0:3]


#squelette = Squelette_3d(f"/data100/fcastillo/RESULT/{snapshots[n]}/{i}_densite_smooth2_c0.1.up.NDskl.S001.a.NDskl")
#xyz = np.array([[p.x, p.y, p.z] for p in squelette.Pointscrit])

taille = []
taille2 = []

Ntot = len(xyz)
print(Ntot)

for b in np.arange(0.5,5.1,0.5) :
    clusters = identify_clusters(xyz, b * (500**3/Ntot)**(1/3))
    taille.append(max(len(cluster) for cluster in clusters))

    clusters_sorted = sorted(clusters, key=len, reverse=True)

    if len(clusters_sorted) >= 2:
        deuxieme = len(clusters_sorted[1])
        taille2.append(deuxieme)

taille = np.array(taille)
taille2 = np.array(taille2)
print(taille, taille/Ntot, taille2, taille2/Ntot)

np.save(f"/data100/fcastillo/RESULT/{snapshots[n]}/{i}_densite_smooth2_c0.1_percolation.txt", taille/Ntot)
np.save(f"/data100/fcastillo/RESULT/{snapshots[n]}/{i}_densite_smooth2_c0.1_percolation_2.txt", taille2/Ntot)