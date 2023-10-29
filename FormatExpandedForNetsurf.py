import pandas as pd
import sys
import os
from os import listdir
from os.path import isfile, join

# Prend tous les fichiers FAA et les réunis dans un même document séparé par motif.
# Permet de fournir un seul document au fichier collab.
# Le résultat s'appelle ConcatenatedFull.faa
if __name__ == '__main__':
    arr = os.listdir("BED_files/Expanded")
    with open("BED_files/Expanded/concatenatedFull.faa", 'w') as fp:
        for file in arr:
            if "Expanded.faa" in file:
                with open("BED_files/Expanded/" + file, 'r') as fd:
                    lines = fd.readlines()
                    motif = file[:(file.find(".faa"))]
                    fp.write("#" + motif + "\n")
                    for number, line in enumerate(lines):
                        fp.write(line)
