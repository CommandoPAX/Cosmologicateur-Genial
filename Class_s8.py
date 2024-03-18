import numpy as np
import matplotlib.pyplot as plt
import yt
import mass_function_library as MFL
from HaloStats import halo_MF
from colossus.cosmology import cosmology
from colossus.lss import mass_function
import os
import json
import re

fichier = open("../../2024/CLASS_NL.dat","r")
k = []
Pk0 = []

for i in range(4):
    fichier.readline()

while 1:
    l = fichier.readline()

    for i in range(10):
        l = l.replace("  "," ")
        l = l.replace("\t", " ")

    if l == "" : break
    l = l.split(" ")
    k.append(float(l[1]))
    Pk0.append(float(l[2].replace("\n","")))

k = np.array(k)
Pk0 = np.array(Pk0)
print(MFL.sigma(k, Pk0, 8.0))
