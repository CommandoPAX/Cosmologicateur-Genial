import os
import sys,getopt

def Recup(argv):
    opts, args = getopt.getopt(argv,"ga")
    
    for opt, arg in opts:
        if opt in ("-a"):
            os.system("scp -r tbruant@obas-hpc.astro.unistra.fr:~/Cosmologicateur-Genial/RESULT/*/*.png ~/RESULT/")
        if opt in ("-o"):
            os.system("scp -r \"tbruant@obas-hpc.astro.unistra.fr:~/Cosmologicateur-Genial/RESULT/* ~/RESULT/")

if __name__=="__main__":
    Recup(sys.argv[1:])