import pandas as pd
import sys
import os
from importlib.metadata import distribution,metadata,version
import biolib
print(biolib.__version__)
import json
from os import listdir
from os.path import isfile, join


# Prend tous les fichiers FAA et les réunis dans un même document séparé par motif.
# Permet de fournir un seul document au fichier collab.
# Le résultat s'appelle concatenatedFull.faa
if __name__ == '__main__':
    arr = os.listdir("BED_files/Expanded")
    with open("BED_files/Expanded/concatenatedFull.faa", 'w') as fp:
        for file in arr:
            if "Expanded.faa" in file and file != "Expanded.faa" and file != "MotifExpanded.faa":
                with open("BED_files/Expanded/" + file, 'r') as fd:
                    lines = fd.readlines()
                    motif = file[:(file.find(".faa"))]
                    fp.write("#" + motif + "\n")
                    for number, line in enumerate(lines):
                        line = line.replace(":", "-")
                        fp.write(line)


    with open("BED_files/Expanded/concatenatedFull.faa", 'r') as fp:
        lines = fp.readlines()
        currentMotif = ""
        i = 0
        os.makedirs('Netsurf_results', exist_ok=True)
        for i in range(len(lines)):
            query = ""
            line = lines[i]
            # line = line.replace(":","=")
            if "#" in line:
                currentMotif = line[1:-1]
                i += 1
                line = lines[i]
                while line[0] != "#" and i < len(lines) - 1:
                    query += line
                    if line[0] != ">":
                        query += "$"
                    i += 1
                    line = lines[i]
                    # line = line.replace(":", "=")

                if "#" not in line:
                    query += line

                querytab = query.split("$")
                print(len(querytab))
                if len(querytab) > 499:
                    i = 1
                    while len(querytab) > 499:
                        part = i
                        querytabPart = querytab[0:250]
                        querytab = querytab[250:]
                        filename = currentMotif + str(part) + ".faa"
                        with open("Netsurf_results/" + filename, 'w') as fw:
                            for j in querytabPart:
                                if len(j) > 0:
                                    # if j[0] == "$":
                                    # j = j[1:]
                                    fw.write(j)  # [:-1])
                        nsp3 = biolib.load('DTU/NetSurfP-3')
                        print(filename)
                        nsp3_results = nsp3.cli(args='-i ' + "Netsurf_results/" + filename)
                        nsp3_results.save_files("Netsurf_results/" + currentMotif + str(part) + "/")
                        i += 1
                else:
                    filename = currentMotif + ".faa"
                    with open("Netsurf_results/" + filename, 'w') as fw:
                        for j in querytab:
                            if len(j) > 0:
                                # if j[0] == "$":
                                # j = j[1:]
                                fw.write(j)  # [:-1])

                    nsp3 = biolib.load('DTU/NetSurfP-3')
                    # Run NetSurfP-3.0 for our query
                    print(currentMotif)
                    print(filename)
                    nsp3_results = nsp3.cli(args='-i ' + "Netsurf_results/" + filename)
                    # Optionally save the results

                nsp3_results.save_files("Netsurf_results/" + currentMotif + "/")

            else:
                i += 1

    print('Traitement des données...')
    for filename in os.listdir('Netsurf_results/'):
        #print(filename)
        #print(os.path.isdir('Netsurf_results/'filename))
        if ".txt" not in filename and ".faa" not in filename and filename != "Expanded" and filename != "MotifExpanded":
            with open('Netsurf_results/' + filename + '/results.json', 'r') as f:
                data = json.load(f)
            # On recupere les données d'interet du fichier Json
            interest = ''
            for e in data:
                structure8 = e['q8']
                structure3 = e['q3']

                for k in range(len(structure8)):
                    # On remplace les élements de structures Coil par un '-' pour plus de
                    # visibilité des hélices
                    if structure8[k] == 'C':
                        structure8 = structure8[:k] + '-' + structure8[k + 1:]
                    if structure3[k] == 'C':
                        structure3 = structure3[:k] + '-' + structure3[k + 1:]

                name = e['desc']
                sequence = e['seq']
                interest += name + '\t'
                interest += structure3 + '\t' + filename + '\n'

            with open('Netsurf_results/' +filename + '.txt', 'w') as w:
                w.write(interest)

    #CutExpanded coordinates
    with open("netsurfHits.txt", 'w') as fp:
        motifList = []
        dictMotifs = {}
        #shortDict = {}
        for file in os.listdir('Netsurf_results/'):
            if "Expanded.txt" in file:
                motif = file[:(file.find("Expanded.txt"))]
                motifList.append(motif)

            elif file[-5].isdigit() and file[-4:] == ".txt":
                motif = file[:(file.find("Expanded"))]
                nextpart = int(file[file.find("Expanded")+8:file.find(".txt")])
                #partList.append(nextpart)
                if dictMotifs.get(motif) is not None:
                    dictMotifs[motif].append(nextpart)
                else:
                    dictMotifs.update({motif: [nextpart]})

        #Concatenate all files
        for motif in dictMotifs.keys():
            with open('Netsurf_results/' + motif + "Expanded.txt", 'w') as motifFile:
                for parts in sorted(dictMotifs[motif]):
                    with open('Netsurf_results/' + motif + "Expanded" + str(parts) + ".txt") as partFile:
                        for partLine in partFile:
                            motifFile.write(partLine)
            motifList.append(motif)

        print(motifList)
        print(dictMotifs)
        seqnameList = []
        motifList2 = []
        structureList = []
        realStartList = []
        realEndList = []
        for motif in motifList:
            fnet = pd.read_csv("Netsurf_results/" +motif + "Expanded.txt", sep="\t",names = ["Name", 'Structure',"Motif"])
            fnet.sort_values(by="Name",inplace=True)

            dfbed = pd.read_csv("BED_files/" + motif + ".bed", sep='\t', names=["Name", "Start", "End","Motif"])
            dfbed.sort_values(by=["Name","Start"],inplace=True)
            i = 0
            j = 0

            for j in range(fnet.shape[0]):
                row = fnet.iloc[j,:]
                line = row["Structure"]
                nameline = row["Name"].split("-")
                name = nameline[0]
                start = nameline[1]
                end = nameline[2]
                seqname = name
                if seqname != dfbed.iloc[i,0]:
                    print("Erreur ici: " + motif + " " + seqname + " " + start+ " " + end)
                    print(dfbed.iloc[i,0])
                    break
                seqnameList.append(seqname)
                realStart = int(dfbed.iloc[i,1])
                realEnd = int(dfbed.iloc[i,2])
                startDiff = realStart - int(start)
                endDiff = int(end) - realEnd
                j = j + 1
                i = i + 1
                motifList2.append(motif)
                structureList.append(line)
                realStartList.append(realStart)
                realEndList.append(realEnd)
                fp.write(seqname + "/" + str(realStart) + "-" + str(realEnd) + "\n")
                fp.write(line + "\n")
                j = j+1

        dict = {"Name":seqnameList,"Start":realStartList,"End": realEndList, "Structure": structureList,"Motif":motifList2}
        dfFull = pd.DataFrame.from_dict(dict)
        dfFull = dfFull.sort_values(by=["Name","Start"])
        dfFull.to_csv("netsurfHelices.txt", sep= ",", index=False, header=False)