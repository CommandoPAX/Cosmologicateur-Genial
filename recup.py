import os
import sys

def Recup(argv):
    opts, args = getopt.getopt(argv,"goa")
    
    for opt, arg in opts:
        if opt in ("-g","-a"):
            POT = True
        elif opt in ("-o"):
            VEL = True
        elif opt in ("-a"):


    os.system()

if __name__=="__main__":
    Recup(sys.argv[1:])