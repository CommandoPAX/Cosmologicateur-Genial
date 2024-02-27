import os
import sys

def Recup(argv):
    opts, args = getopt.getopt(argv,"goa")
    
    for opt, arg in opts:
        if opt in ("-g","-a"):
            os.system()
        if opt in ("-o","-a"):
            os.system()


    os.system()

if __name__=="__main__":
    Recup(sys.argv[1:])