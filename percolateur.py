from skeletonnateur_hpc import*
import sys
import numpy as np

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

n = 0 #int(sys.argv[1])
i = 4 #int(sys.argv[2])

snapshots = ["benchM","NG_F500","G_m500","NG_F500_m500","NG_Fminus500","NG_Fminus500_m500", "G_ViVi", "NG_ViVi", "NG_Fminus500_ViVi"]
#labels = ["LCDM", "fnl = -500", "m = 500 eV", "WDM & fnl = -500", "fnl = 500", "WDM & fnl = 500"]

#snapshots = ["NEDE","NsPNG_F500","NsPNG_F1000","NsPNG_F1833","NsPNG_EDE_F500","NsPNG_EDE_F1000","NsPNG_EDE_F1833"]

z= [15,12, 10, 8, 5,3,1,0.5,0.25,0]
indices_z = [5,6,8,9]



squelette = Squelette_3d(f"/data100/fcastillo/RESULT/{snapshots[n]}/{i}_densite_smooth2_c0.1.up.NDskl.S001.a.NDskl")
xyz = np.array([[p.x, p.y, p.z] for p in squelette.Pointscrit])
taille = []

for b in np.arange(0.2,1.2,0.05) :
    clusters = identify_clusters(xyz, b * (500**3/len(squelette.PointsCrit))**(1/3))

    taille.append(max(len(cluster) for cluster in clusters))

plt.plot(np.arange(0.2,1.2,0.05), np.array(taille)/len(squelette.PointsCrit))
plt.savefig("percolation.pdf")