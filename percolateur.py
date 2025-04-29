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

n = int(sys.argv[1])
i = int(sys.argv[2])

snapshots = ["NEDE","NsPNG_F500","NsPNG_F1000","NsPNG_F1833","NsPNG_EDE_F500","NsPNG_EDE_F1000","NsPNG_EDE_F1833"]

z= [15,12, 10, 8, 5,3,1,0.5,0.25,0]
indices_z = [5,6,8,9]



squelette = Squelette_3d(f"/data100/fcastillo/RESULT/{snapshots[n]}/{i}_densite_smooth2_c0.1.up.NDskl.S001.a.NDskl")
xyz = np.array([[p.x, p.y, p.z] for p in squelette.Pointscrit])
taille = []
taille2 = []

for b in np.arange(0.6,1.4,0.01) :
    clusters = identify_clusters(xyz, b * (500**3/len(squelette.Pointscrit))**(1/3))
    taille.append(max(len(cluster) for cluster in clusters))

    clusters_sorted = sorted(clusters, key=len, reverse=True)

    if len(clusters_sorted) >= 2:
        deuxieme = len(clusters_sorted[1])
        taille2.append(len(deuxieme))

np.save(f"/data100/fcastillo/RESULT/{snapshots[n]}/{i}_densite_smooth2_c0.1_percolation.txt", np.array(taille)/len(squelette.Pointscrit))
np.save(f"/data100/fcastillo/RESULT/{snapshots[n]}/{i}_densite_smooth2_c0.1_percolation_2.txt", np.array(taille2)/len(squelette.Pointscrit))