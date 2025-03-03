import numpy as np
from math import sqrt
import sys

class Simplexe ():
    def __init__ (self, ndim, points) :
        self.ndim = ndim
        self.points = points

    def Voulme (self) :
        pass

    def Surface (self):
        if self.ndim == 2:
            AB = np.array([self.points[1].x - self.points[0].x,
                           self.points[1].y - self.points[0].y,
                           self.points[1].z - self.points[0].z] )
            AC = np.array([self.points[2].x - self.points[0].x,
                           self.points[2].y - self.points[0].y,
                           self.points[2].z - self.points[0].z] )
            
            vect = np.array([AB[1]*AC[2] - AB[2]*AC[1],
                             AB[2]*AC[0] - AC[2]*AB[0],
                             AB[0]*AC[1] - AB[1]*AC[0]])
            
            return sqrt(vect[0]**2 + vect[1]**2 + vect[2] **2)

class Vertex ():
    def __init__(self, x, y, z):
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)

class Net ():
    def __init__ (self, path) :
        self.path = path
        self.fichier = open(path, "r")
        fichier = self.fichier

        fichier.readline()
        self.ndim = int(fichier.readline()[:-1])
        fichier.readline()
        fichier.readline()
        self.num_vertices = int(fichier.readline()[:-1])

        self.list_vertices = []
        for i in range(self.num_vertices):
            coords = (fichier.readline()[:-1]).split(" ")
            self.list_vertices.append(Vertex(coords[0], coords[1], coords[2]))

        self.liste_simplexes = []

        while True :
            ligne = fichier.readline()
            if "[" in ligne : break
            data_simplex = (ligne[:-1]).split(" ")
            T = data_simplex[0]
            N = data_simplex[1]

            for i in range(N):
                points = (fichier.readline()[:-1]).split(" ")
                self.liste_simplexes.append(Simplexe(int(T),points))

        fichier.close()

if __name__ == "__main__" :
    n = int(sys.argv[1])
    i = int(sys.argv[2])

    snapshots = ["benchM","NG_F500","G_m500","NG_F500_m500","NG_Fminus500","NG_Fminus500_m500"]
    labels = ["LCDM", "fnl = -500", "m = 500 eV", "WDM & fnl = -500", "fnl = 500", "WDM & fnl = 500"]


    reseau = Net(f"/data100/fcastillo/RESULT/{snapshots[n]}/{i}_densite_0_c0.1_manifolds_J01a.NDnet.S001.a.NDnet")

    surfaces = []
    for triangle in reseau.liste_simplexes :
        surfaces.append(triangle.Surface())

    np.save(f"/data100/fcastillo/RESULT/{snapshots[n]}/{i}_densite_0_c0.1_surf_murs.txt", np.array(surfaces))


