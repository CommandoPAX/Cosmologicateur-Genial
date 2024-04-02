import os

for root, dirs, files in os.walk("./RESULT/"):

    for dir in dirs :

        for root2,dirs2,files2 in os.walk("./RESULT/"+dir):

            max = 0

            for i in range(7):
                for file in files2 :
                    if str(i)+"_" in file and i > max :
                        max = i

            if max == 6 :
                for file in files2 :
                    if "2_" in file or "3_" in file or "4_" in file :
                        os.system("rm \"./RESULT/"+dir+"/"+file+"\"")

                for file in files2 :
                    path = "RESULT/"+dir+"/"

                    if "6_" in file :
                        nouveau = "3_"+file[2:]
                        os.system("mv \""+path+file+ "\" \""+path+nouveau+"\"")
                    if "5_" in file :
                        nouveau = "2_"+file[2:]
                        os.system("mv "+file+ " "+nouveau)

            if max == 5 :
                for file in files2 :
                    if "2_" in file or "3_" in file:
                        os.system("rm \"./RESULT/"+dir+"/"+file+"\"")

                for file in files2 :
                    path = "RESULT/"+dir+"/"

                    if "5_" in file :
                        nouveau = "3_"+file[2:]
                        os.system("mv \""+path+file+ "\" \""+path+nouveau+"\"")
                    if "4_" in file :
                        nouveau = "2_"+file[2:]
                        os.system("mv "+file+ " "+nouveau)

            if max == 4 :
                for file in files2 :
                    if "2_" in file :
                        os.system("rm \"./RESULT/"+dir+"/"+file+"\"")

                for file in files2 :
                    path = "RESULT/"+dir+"/"

                    if "4_" in file :
                        nouveau = "3_"+file[2:]
                        os.system("mv \""+path+file+ "\" \""+path+nouveau+"\"")
                    if "3_" in file :
                        nouveau = "2_"+file[2:]
                        os.system("mv "+file+ " "+nouveau)

for root, dirs, files in os.walk("./RESULT/"):

    for dir in dirs :

        for root2,dirs2,files in os.walk("./RESULT/"+dir):

            for dir2 in dirs2 :
                print (dir+" -> "+dir2)