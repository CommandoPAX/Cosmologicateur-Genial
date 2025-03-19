import os
from pdf2image import convert_from_path
from wand.image import Image

for root, dirs, files in os.walk("./WDM"):
    for file in files:
        with Image(filename=root+"/"+file) as img:
            if ".pdf" in file : img.save(filename="./PNG/"+file[:-4]+".png")
