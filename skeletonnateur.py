import numpy as np
import os
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from astropy.io import fits

class PointCrit ():
    def __init__ (self, type, x, y, z, value, pairid, boundary, nfil, data_fil):

        self.type = type
        self.x = x
        self.y = y
        self.z = z
        self.value = value
        self.pairid = pairid
        self.boundary = boundary
        self.nfil = nfil
        self.data_fil = data_fil

class Filament ():
    def __init__ (self, cp1, cp2, nsamp, data_samp):
        self.cp1 = cp1
        self.cp2 = cp2
        self.nsamp = nsamp
        self.data_samp = data_samp

class Squelette():
    def __init__ (self, path):

        self.path = path
        fichier = open(self.path,"r")

        fichier.readline()
        self.ndim = fichier.readline()
        self.comment = fichier.readline()
        self.bbox = fichier.readline()

        fichier.close()

        self.Critical_points()

    def Critical_points (self) :

        self.fichier = open(self.path,"r")
        fichier = self.fichier

        while 1 :
            ligne = fichier.readline()
            if "[CRITICAL POINTS]" in ligne : break

        self.ncrit = int(fichier.readline()[:-1])
        self.Pointscrit = []
        for nc in range(self.ncrit) :
            info = fichier.readline()[:-1]
            infos = info.split(" ")
            type_ = infos[0]
            x = infos[1]
            y = infos[2]
            z = infos[3]
            value = infos[4]
            pairid = infos[5]
            boundary = infos[6]
            nfil = int(fichier.readline()[1:-1])

            data_fil = []

            for n in range(nfil) :
                info_fil = fichier.readline()[1:-1]
                infos_fil = info_fil.split(" ")
                data_fil.append(infos_fil)


            self.Pointscrit.append(PointCrit(type_, x, y, z, value, pairid, boundary, nfil, data_fil))
     
        if "[FILAMENTS]" in fichier.readline() :
            self.filaments()

    def filaments(self):
        fichier = self.fichier

        self.nfil = int(fichier.readline()[:-1])
        self.liste_filaments = []

        for nf in range(self.nfil):
            infos = (fichier.readline()[:-1]).split(" ")
            cp1 = infos[0]
            cp2 = infos[1]
            nsamp = int(infos[2])
            data_samp = []
            for ns in range(nsamp):
                info_samp = (fichier.readline()[1:-1]).split(" ")
                data_samp.append(info_samp)

            self.liste_filaments.append(Filament(cp1, cp2, nsamp, data_samp))

        if "[CRITICAL POINTS DATA]" in fichier.readline():
            self.Critical_points_data()

    def Critical_points_data(self):
        fichier = self.fichier

        self.num_fields = int(fichier.readline()[:-1])
        self.nom_fields = []
        for nom in range(self.num_fields):
            self.nom_fields.append(fichier.readline()[:-1])

        for pc in self.Pointscrit :
            pc.value_fields = (fichier.readline()[:-1]).split(" ")

        if "[FILAMENTS DATA]" in fichier.readline() :
            self.Filaments_data()

    def Filaments_data (self):
        fichier = self.fichier

        self.num_fields_fil = int(fichier.readline()[:-1])
        self.nom_fields_fil = []
        for nom in range(self.num_fields_fil):
            self.nom_fields_fil.append(fichier.readline()[:-1])
            
        for fil in self.liste_filaments :
            fil.value_fields = (fichier.readline()[:-1]).split(" ")

        fichier.close()

def plot_filaments(squelette):
    ax = plt.figure().add_subplot(projection='3d')
    i = 0
    j = 0
    color = ["red", "blue", "green", "orange", "black"]
    for fil in squelette.liste_filaments :
        j+=1
        for n in range(len(fil.data_samp)-1):
            p0 = fil.data_samp[n]
            p1 = fil.data_samp[n+1]

            x0 = float(p0[0])
            y0 = float(p0[1])
            z0 = float(p0[2])

            x1 = float(p1[0])
            y1 = float(p1[1])
            z1 = float(p1[2])

            d = (x0 - x1)**2 + (y0 - y1)**2 + (z0-z1)**2

            if ((z0>=0 and z0 <=4) or (z1>=0 and z1 <=4)) and ((x0>=0 and x0 <=4) or (x1>=0 and x1 <=4)) and ((y0>=0 and y0 <=2.5) or (y1>=0 and y1 <=4)) and d <25:
                i+=1

                ax.plot([x0,x1],[y0,y1],[z0,z1],color="blue")
                ax.scatter([x0,x1],[y0,y1],[z0,z1],color="blue")

        if i == 100 : break

    ax.set_zlim(0,4)
    ax.set_xlim(0,4)
    ax.set_ylim(0,4)
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")

    plt.show()

def plot_filaments_2d(squelette, data):

    #fits_image_filename = fits.util.get_testdata_filepath(data)

    hdul = fits.open(data)
    delta = hdul[0].data
    print(np.shape(delta))

    plt.imshow(delta[:, :, 0] + delta[:, :, 1] + delta[:, :, 2] + delta[:, :, 3])

    i = 0
    j = 0
    color = ["red", "blue", "green", "orange", "black"]
    for fil in squelette.liste_filaments :
        j+=1
        for n in range(len(fil.data_samp)-1):
            p0 = fil.data_samp[n]
            p1 = fil.data_samp[n+1]

            x0 = float(p0[0])
            y0 = float(p0[1])
            z0 = float(p0[2])

            x1 = float(p1[0])
            y1 = float(p1[1])
            z1 = float(p1[2])

            d = (x0 - x1)**2 + (y0 - y1)**2 + (z0-z1)**2

            if ((z0>=0 and z0 <=4) or (z1>=0 and z1 <=4)) and d <2000:
                i+=1

                plt.plot([x0,x1],[y0,y1],color="blue")
                #plt.scatter([x0,x1],[y0,y1],color="blue")

        #if i >= 5000 : break



    plt.show()


if __name__ == "__main__" :
    squelette = Squelette("../gaussien_densite.fits_c1.up.NDskl.S001.a.NDskl")
    data = "../gaussien_densite.fits"
    plot_filaments_2d(squelette, data)
