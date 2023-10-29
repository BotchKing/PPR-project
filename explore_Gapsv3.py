import sys
import os
from os import listdir
from os.path import isfile, join
import pandas as pd
import numpy

# On veut examiner nos séquences de confiance et trouver nos candidats dans les "gaps"
# Une fois les gaps identifiés, on veut analyser leur structure
# prendre nos fichiers BED révisés avec seulement les séquences Trusted.

def checkGaps(df):

    f_length = df.shape[0]
    i = 0
    df = df.drop([0])
    while (i < f_length-1):
        row = df.iloc[i,:]
        name1 = row["Name"]
        evalue = row["E_value"]
        cevalue = row["c_Evalue"]
        ievalue = row["i_Evalue"]
        end1 = int(row["End"])
        tlen = int(row["tlen"])
        j= i+1
        row2 = df.iloc[j,:]
        name2 = row2["Name"]
        start2 = int(row2["Start"])
        score = ""
        
        if (name1 == name2):
            if (((start2-end1)<=40) and ((start2-end1)>=25)):
               newrow = [name1,tlen,(end1+1),(start2-1),'Candidate++',score]

            else:
                i = i+1
                continue
            newrowdf = pd.DataFrame([newrow],columns=['Name',"tlen",'Start','End','Motif','score2'])
            df2 =pd.concat([df,newrowdf],axis=0,join='outer',ignore_index=True)
            df= df2

        i = i+1
    return df


if __name__ == '__main__':
    identifiers = []
    filename = sys.argv[1]
    df = pd.read_csv(filename, sep="\t",names=["Name", "tlen", "Motif", "qlen", "E_value",
               "c_Evalue", "i_Evalue", "Start", "End","score2"],index_col=False)
    #while (df.shape[1] >6):
    #    df.pop(df.columns[-1])

    df3 = checkGaps(df)
    df3 = df3.sort_values(by=['Name','Start','End'])
    df3.to_csv(str(filename)[:-4]+"_AndCandidates.bed", sep='\t',index=None)
    print('added the candidate motifs in 25-40 gaps')

