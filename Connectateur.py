# coding: utf-8

import socket

hote = socket.gethostbyname(socket.gethostname())

port = 15555

socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
socket.connect((hote, port))
print ("Connection on {}".format(port))

socket.sendto("Bonjour de l'HPC".encode(),(hote,port))

print ("Close")
socket.close()
