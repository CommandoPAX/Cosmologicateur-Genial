import os
import getpass
import subprocess
import pexpect

def Recup(argv):
    opts, args = getopt.getopt(argv,"gac")
    
    for opt, arg in opts:
        if opt in ("-g"):
            commande = "scp -r tbruant@obas-hpc.astro.unistra.fr:~/Cosmologicateur-Genial/RESULT/*/*.png ~/RESULT/"
        if opt in ("-a"):
            commande = "scp -r tbruant@obas-hpc.astro.unistra.fr:~/Cosmologicateur-Genial/RESULT/* ~/RESULT/"
        
        if opt in ("-c"):
            commande = "scp -r tbruant@obas-hpc.astro.unistra.fr:~/Cosmologicateur-Genial/*.py ~/Cosmologicateur-genial/"
        
        var_password = str("gaeK6Oafai1eN6e")

        var_child = pexpect.spawn(commande)
        i = var_child.expect(["password:", pexpect.EOF])

        if i==0: # send password            
            var_child.sendline(var_password)
            var_child.expect(pexpect.EOF)
        elif i==1: 
            print("Got the key or connection timeout")
            pass

if __name__=="__main__":
    Recup(sys.argv[1:])