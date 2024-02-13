# coding: utf-8
##### Test connection à un serveur local pour transmettre les données

import socket

hote = "localhost"
port = 15555

socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
socket.connect((hote, port))
print ("Connection on {}".format(port))

socket.sendto("Bonjour de l'HPC".encode(),(hote,port))

print ("Close")
socket.close()
