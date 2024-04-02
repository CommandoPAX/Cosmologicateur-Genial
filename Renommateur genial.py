import os

for root, dirs, files in os.walk("./RESULT/"):

    for dir in dirs :

        for root2,dirs2,files in os.walk("./RESULT/"+dir):

            if "output_00005" in dirs2 :
                os.system("rm -r \"./RESULT/"+ dir+"/output_00002\"")
                os.system("mv \"./RESULT/"+ dir+"/output_00003\" \"./RESULT/"+ dir+"/output_00002\"" )
                os.system("mv \"./RESULT/"+ dir+"/output_00005\" \"./RESULT/"+ dir+"/output_00003\"" )
           
            
for root, dirs, files in os.walk("./RESULT/"):

    for dir in dirs :

        for root2,dirs2,files in os.walk("./RESULT/"+dir):

            for dir2 in dirs2 :
                print (dir+" -> "+dir2)