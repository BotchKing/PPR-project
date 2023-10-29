import sys
import os
import math
import pandas as pd
import numpy as np

#Il faut appeler ce script du dossier qui contient tous les dossiers de résultats d'alphafold
#Appelé avec une liste de séquences avec tous leurs hits. Va chercher les fichiers qui contiennent les hélices
#produits par ChimeraXScript.py. On veut extraire les patterns d'hélices dans les régions hits

def removeWrongStructures(originalDf,dataframeAlpha,dataframeNetsurf):
    notMotifsIndex = []
    # Time to decide our conditions for accepting a candidate as a helix-turn-helix motif
    # For now, let's say a candidate is accepted if it does not feature a "-1" on either h1,turn or h2 feature
    os.makedirs("NotMotifs", exist_ok=True)
    structMerged = pd.merge(dataframeAlpha,dataframeNetsurf, on=["Name","Sequence Start", "Sequence End"], how="right")
    structMerged["Final Structure"] = structMerged.apply(lambda row:row["Structure_y"] if pd.isna(row["Structure_x"]) else row["Structure_x"],axis=1)
    structMerged["#H1"] = structMerged.apply(
        lambda row: row["H1 Count_y"] if pd.isna(row["H1 Count_x"]) else row["H1 Count_x"], axis=1)
    structMerged["#H2"] = structMerged.apply(
        lambda row: row["H2 Count_y"] if pd.isna(row["H2 Count_x"]) else row["H2 Count_x"], axis=1)
    structMerged["Motif_xy"] = structMerged.apply(
        lambda row: row["Motif_y"] if pd.isna(row["Motif_x"]) else row["Motif_x"], axis=1)

    structMerged = structMerged[["Name","Sequence Start","Sequence End", "Final Structure", "#H1", "#H2","Motif_xy"]]
    structMerged["H1Criteria"] = structMerged["#H1"].lt(7)
    structMerged["H2Criteria"] = structMerged["#H2"].lt(7)
    keeplist = []
    for row in range(structMerged.shape[0]):
        keepvalue = structMerged.iloc[row,-2] or structMerged.iloc[row,-1]
        keeplist.append(keepvalue)
    structMerged["Keep"] = keeplist
    #structMerged["Keep"] = pd.np.where(structMerged[["#H1"]] < 7,"No", pd.np.where(structMerged[["#H2"]] < 7,"No","Yes"))
    structMerged.to_csv("NotMotifs/FullStructures.bed", sep="\t", index=False)

    droppedDf = structMerged.loc[(structMerged["#H1"] < 7) | (structMerged["#H2"] < 7)]
    
    droppedDf.to_csv("NotMotifs/AllRemoved.bed", sep="\t")
    buildDf = droppedDf[["Name", "Sequence Start", "Sequence End", "Motif_xy","Keep"]]

    buildDf = buildDf.rename(columns={"Sequence Start":"Start","Sequence End":"End", "Motif_xy": "Motif"})
    buildDf.to_csv("NotMotifs/FullDrop.bed", sep="\t",index=False)

    endDf = structMerged[structMerged["Keep"] != True]
    endDf.to_csv("NotMotifs/MergedCleaned.bed", sep="\t", index=False)
    dfOut = endDf[["Name", "Sequence Start", "Sequence End", "Motif_xy",]]
    dfOut = dfOut.rename(columns={"Sequence Start": "Start", "Sequence End": "End", "Motif_xy": "Motif"})
    dfOut.to_csv("NotMotifs/CleanedList.bed", sep="\t",index=False)

def composition(df, filename):
    with open(filename, 'w') as fw:
        listHits = []
        for i in range(df.shape[0]):
            seqName = df.iloc[i, 0]
            seqStart = df.iloc[i, 1]
            seqEnd = df.iloc[i, 2]
            linet = str(df.iloc[i, 3])
            motif = str(df.iloc[i, 4])
            fw.write(seqName + "\n")
            fw.write(linet + "\t" + motif + "\n")
            index = 0
            inH1 = False
            inTurn = False
            inH2 = False
            h1count = 0
            h1end = -1
            turnstart = -1
            turnend = -1
            h2count = 0
            h2end = len(linet)
            previous= -1
            secondprevious = -1
            while index < (len(linet) - 2):
                letter = linet[index]
                if letter == "-":
                    letter = "+"
                if index <= 10:
                    if letter == "H":
                        h1count += 1
                        if secondprevious != -1:
                            secondprevious = previous
                        previous = "H"
                        index += 1

                    else:
                        if secondprevious != -1:
                            secondprevious = previous
                        previous = "+"
                        index += 1
                else:
                    #Above 13 index, we start looking at the previous letters
                    if (secondprevious == "+") and (previous == "+") and (letter == "+") and (not inTurn):
                        #We are in the turn, if the number of helices in the first helix are above 10, we keep it
                        inTurn = True

                    elif inTurn:
                        if letter == "H":
                            h2count += 1
                            index += 1
                        else:
                            index += 1
                    else:
                        if letter == "H":
                            h1count += 1
                            secondprevious = previous
                            previous = letter
                            index += 1
                        else:
                            secondprevious = previous
                            previous = letter
                            index += 1
            listHits.append([seqName, seqStart, seqEnd, linet, h1count, h2count, motif])
            fw.write("h1:" + str(h1count) + "\n")
            fw.write("h2: " + str(h2count) + "\n")
            fw.write("............................................ \n")


    dfOut2 = pd.DataFrame(listHits,
                          columns=['Name','Sequence Start','Sequence End','Structure', 'H1 Count', 'H2 Count', "Motif"])

    return dfOut2




if __name__ == '__main__':
    #Give the file containing the motifs start and end. Bed format
    file = sys.argv[1]
    #Give the file containing the netsurf results (NetsurfHelices.txt)
    netsurf = sys.argv[2]
    arr = os.listdir("AlphafoldData")
    df = pd.read_csv(file,sep='\t',names=["Name","Start","End","Motif"],header=0)
    netsurfDf = pd.read_csv(netsurf, sep=',', names=["Name", "Start","End","Structure", "Motif"])
    i = 0
    intervalList = []
    listdict = {}
    namesList= []
    startList=[]
    endList=[]
    motifList=[]
    structureList = []
    ngroupName = df.groupby(["Name"]).ngroups
    with open("VisualHelices.txt", 'w') as fp:
        for group in df.groupby(["Name"]):
            for directory in arr:

                if(os.path.isdir("AlphafoldData/" + str(directory))):
                    datapath = "AlphafoldData" + str(directory)
                    if (os.path.exists(datapath) and datapath == group[0]):
                        with open(str(directory)+ "/" + str(directory) + ".txt",'r') as fd:
                            name = group[0]
                            fp.write(name + "\t" + "\n")

                            lines = fd.readlines()

                            for i in range(group[1].shape[0]):
                                j = 0
                                structure = ""
                                hitstart = group[1].iloc[i, 1]
                                hitend = group[1].iloc[i, 2]
                                DFinterval = pd.Interval(left=int(hitstart),
                                                         right=int(hitend), closed="both")
                                for number, line in enumerate(lines):
                                    if ";" in line:
                                        numbers = line.split(";")
                                        hstart = int(numbers[0])
                                        hend = int(numbers[1])
                                        Hinterval = pd.Interval(left=hstart, right=hend, closed="both")

                                        #if DFinterval.overlaps(Hinterval):
                                            #listOut.append(hinterval)

                                        # we show viually where the helices are in a new file
                                        while j < hend:
                                            isInHit = j in DFinterval

                                            # "-" when not helix AND not in a reported hit
                                            if (j < hstart) and (isInHit == False):
                                                fp.write("-")
                                                j += 1

                                            # "+" when not helix AND  in a reported hit
                                            elif (j < hstart) and isInHit:
                                                fp.write("+")
                                                structure += "+"
                                                j += 1

                                            # "H" when in helix AND  in a reported hit
                                            elif (j >= hstart) and isInHit:
                                                fp.write("H")
                                                structure += "H"
                                                j += 1

                                            # "h" when in helix AND  in not a reported hit
                                            # elif (j >= hstart) and (isInHit == False):
                                            else:
                                                fp.write("h")
                                                j += 1

                                while j < hitend:
                                    structure += "+"
                                    j += 1
                                motif = group[1].iloc[i, 3]
                                namesList.append(name)
                                startList.append(hitstart)
                                endList.append(hitend)
                                structureList.append(structure)
                                motifList.append(motif)
                        fp.write("\n")


    StructDict = {"Name":namesList,"Start":startList,"End":endList,"Structure":structureList,"Motif": motifList}
    hitdf = pd.DataFrame.from_dict(StructDict)
    hitdf.to_csv("AlphafoldHitHelices.txt", sep=",", index=False, header=False)


    dfAlpha = composition(hitdf,"CompositionAlphafold.txt")
    dfAlpha.to_csv("AlphafoldAllData.csv",index=False)

    dfNet = composition(netsurfDf,"CompositionNetsurf.txt")
    removeWrongStructures(df, dfAlpha, dfNet)

