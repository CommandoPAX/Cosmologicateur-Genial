import os

for root, dirs, files in os.walk("./RESULT/"):

    for dir in dirs :

        for root2,dirs2,files in os.walk("./RESULT/"+dir):
            for dir2 in dirs2 :
                print (dir+" -> "+dir2)