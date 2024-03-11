import os
import sys,getopt
from subprocess import Popen, PIPE
import time

"""p = Popen(['sudo', '-S', 'ls'], stdin=PIPE, stderr=PIPE, stdout=PIPE, text=True)
prompt = p.communicate("password" + '\n')
output = prompt[0]"""

#print(output)

def Recup(argv):
    opts, args = getopt.getopt(argv,"gac")
    
    for opt, arg in opts:
        if opt in ("-g"):
            p = Popen(['scp', '-r', '\"tbruant@obas-hpc.astro.unistra.fr:~/Cosmologicateur-Genial/RESULT/*/*.png\"','~/RESULT/'], stdin=PIPE, stderr=PIPE, stdout=PIPE, text=True)
            
            time.sleep(2)

            p.stdin.write("gaeK6Oafai1eN6e" + '\n')
            p.stdin.flush()
            #os.system("scp -r tbruant@obas-hpc.astro.unistra.fr:~/Cosmologicateur-Genial/RESULT/*/*.png ~/RESULT/")
        if opt in ("-a"):
            #os.system("scp -r \"tbruant@obas-hpc.astro.unistra.fr:~/Cosmologicateur-Genial/RESULT/*\" ~/RESULT/")
        
            p = Popen(['scp', '-r', '\"tbruant@obas-hpc.astro.unistra.fr:~/Cosmologicateur-Genial/RESULT/*\"','~/RESULT/'], stdin=PIPE, stderr=PIPE, stdout=PIPE, text=True)
            
            time.sleep(2)

            
            p.stdin.write("gaeK6Oafai1eN6e" + '\n')
            p.stdin.flush()
        
        if opt in ("-c"):
            p = Popen(['scp', '-r', '\"tbruant@obas-hpc.astro.unistra.fr:~/Cosmologicateur-Genial/*.py\"','~/Cosmologicateur-genial/'], stdin=PIPE, stderr=PIPE, stdout=PIPE, text=True)
            
            time.sleep(2)
            
            p.stdin.write("gaeK6Oafai1eN6e" + '\n')
            p.stdin.flush()
            #os.system("scp -r \"tbruant@obas-hpc.astro.unistra.fr:~/Cosmologicateur-Genial/*.py\" ~/Cosmologicateur-genial/")
        


if __name__=="__main__":
    Recup(sys.argv[1:])