import numpy as np
import os
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from astropy.io import fits
from math import*

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

    def len (self) :

        self.l = 0

        for n in range(len(self.data_samp)-1):
            p0 = self.data_samp[n]
            p1 = self.data_samp[n+1]

            for ndim in range(len(p0)):
                self.l += sqrt((float(p0[ndim]) - float(p1[ndim]))**2) * 500/512


class Squelette_2d():
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
            value = infos[3]
            pairid = infos[4]
            boundary = infos[5]
            nfil = int(fichier.readline()[1:-1])

            data_fil = []

            for n in range(nfil) :
                info_fil = fichier.readline()[1:-1]
                infos_fil = info_fil.split(" ")
                data_fil.append(infos_fil)


            self.Pointscrit.append(PointCrit(type_, x, y, 0, value, pairid, boundary, nfil, data_fil))
     
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

    def plot(self, axe, field) :

        if True :
            im = axe.imshow(field, origin="lower", vmax = 10, vmin = -2)
            plt.colorbar(im)

        axe.set_xlim(300,400)
        axe.set_ylim(300,400)


        for fil in self.liste_filaments :
        
            for n in range(len(fil.data_samp)-1):
                p0 = fil.data_samp[n]
                p1 = fil.data_samp[n+1]

                x0 = float(p0[0])
                y0 = float(p0[1])

                x1 = float(p1[0])
                y1 = float(p1[1])

                d = (x0 - x1)**2 + (y0 - y1)**2

                if d <2000 :

                    axe.plot([x0,x1],[y0,y1],color="white", linewidth=1)

 

class Squelette_3d():
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

            if d<25 : #((z0>=0 and z0 <=4) or (z1>=0 and z1 <=4)) and ((x0>=0 and x0 <=4) or (x1>=0 and x1 <=4)) and ((y0>=0 and y0 <=2.5) or (y1>=0 and y1 <=4)) and d <25:
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


    plt.imshow((np.sum(delta[:,:,0:2], axis=2)), origin="lower", vmax = 10)
    plt.colorbar()

    data_slice = delta[:,:,0:2]

    """header = fits.Header()
    header['COMMENT'] = 'Slice'
    header['NAXIS'] = 3  
    header['NAXIS1'] = delta.shape[0]  
    header['NAXIS2'] = delta.shape[1]  
    header['NAXIS3'] = delta.shape[2] 

    hdu = fits.PrimaryHDU(data=data_slice, header=header)

    hdu.writeto(f"slice_densite.fits", overwrite=True)"""



    i = 0
    j = 0
    color = ["red", "blue", "green", "orange", "black"]
    for fil in squelette.liste_filaments :
        j+=1
        if True:#float(fil.value_fields[0])> 1:
            for n in range(len(fil.data_samp)-1):
                p0 = fil.data_samp[n]
                p1 = fil.data_samp[n+1]

                z0 = float(p0[0])
                x0 = float(p0[1])
                y0 = float(p0[2])

                z1 = float(p1[0])
                x1 = float(p1[1])
                y1 = float(p1[2])

                d = (x0 - x1)**2 + (y0 - y1)**2 + (z0-z1)**2

                if d <2000 and z1 <3 and z0 <3:
                    i+=1

                    plt.plot([x0,x1],[y0,y1],color="white", linewidth=1)
                    #plt.scatter([x0,x1],[y0,y1],color="blue")

        #if i >= 5000 : break



    plt.show()

def PDF_len_filaments (axes, squelette, couleur="blue", ls="-") :
    longueurs = []
    for fil in squelette.liste_filaments :
        fil.len()
        longueurs.append(fil.l)

    hist = np.histogram(np.array(longueurs), density= True, range = [0, 10], bins=10)
    axes.plot(hist [0], color= couleur, ls = ls)
    axes.set_xlabel("longueur [Mpc / h]")
    axes.set_ylabel("Probabilite")
    axes.set_ylim(0,0.3)




if __name__ == "__main__" :

    Redshifts = {
        0 : 32,
        1 : 3,
        2 : 1,
        3 : 0.25,
        4 : 0
    }
    
    snapshots = ["benchM","NG_F500","G_m500","NG_F500_m500","NG_Fminus500","NG_Fminus500_m500"]

    ls = ["-", "-", "-.", "--", "-", "--"]
    couleurs = ["blue", "orange", "green", "fuchsia", "orange", "fuchsia"]

    for i in range(1,5):
        plt.subplot(2,2,i)

        axes = plt.gca()

        axes.title.set_text (f"z = {Redshifts[i]}")

        for j in range(6):


            squelette = Squelette_2d(f"/data100/fcastillo/RESULT/{snapshots[j]}/{i}_densite_slice.fits_c1.up.NDskl.S001.a.NDskl")
            data = f"../slice/{i}_densite_slice.fits"
            #field_fichier = fits.open(data)
            #field = field_fichier[0].data
            #squelette.plot(axes, field = field)
            PDF_len_filaments(axes, squelette, couleur = couleurs[j], ls = ls[j])
            #field_fichier.close()

    plt.savefig("len.png")

