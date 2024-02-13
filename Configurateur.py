##### Le configurateur génial
##### Automatise la configuration de Ramses et monofonic

import os

monofonic= open("../monofonic/monofonic.conf","r")
ramses = open("../ramses/namelist/ramses.nml","r")

ngrid = input("Nombre de cellules : 2^")
print("n^"+ngrid +" = "+str(2**int(ngrid)))
taille = input("Taille de la boite : ")

output_monofonic = ""
while 1 :
    ligne = monofonic.readline()
    if ligne == "" : break
    for i in range(10): ligne = ligne.replace("  "," ")
    ligne = ligne.split(" ")
    if ligne[0] == "GridRes" : ligne[2] = str(2**int(ngrid))
    if ligne[0] == "BoxLength" : ligne[2] = taille

    i  = 0
    for k in ligne :
        output_monofonic += k
        if not "\n" in k and not i== 0: 
            output_monofonic+=" "
        i = 1
monofonic.close()


output_ramses = ""
while 1 :
    ligne = ramses.readline()
    if ligne == "" : break
    for i in range(10): ligne = ligne.replace("  "," ")
    ligne = ligne.split("=")
    if ligne[0] == "levelmin" : ligne[1] = ngrid+"\n"
    if ligne[0] == "levelmax" : ligne[1] = str(int(ngrid)+6)+"\n"

    for k in ligne :
        output_ramses += k
        if not "\n" in k : 
            output_ramses+="="
ramses.close()

output_monofonic = output_monofonic.replace("="," = ")
output_monofonic = output_monofonic.replace("#"," # ")
output_monofonic = output_monofonic.replace("# #"," # ")

autre = input("Changer d'autres paramètres ? (o/n)")

monofonic = open("../monofonic/monofonic.conf","w")
monofonic.write(output_monofonic)
monofonic.close()

ramses = open("../ramses/namelist/ramses.nml","w")
ramses.write(output_ramses)
ramses.close()

if autre == "o" :
    os.system("nano ../monofonic/monofonic.conf")
    os.system("nano ../ramses/namelist/ramses.nml")

os.system("../monofonic/build/monofonIC ../monofonic/monofonic.conf")
os.system("cp ../monofonic/ics_ramses/ic_poscx ../monofonic/ics_ramses/ic_deltab")
os.system("cd ../ramses/bin; make clean")
os.system("cd ../ramses/bin; make NDIM=3 ")
os.system("../ramses/bin/ramses3d ../ramses/namelist/ramses.nml")