import os
from os import listdir
from os.path import isfile, join
import pandas as pd
import numpy as np
import sys
#Permet de comparer deux tables de listes PPR, fournit par la fonction extract_motifs.py.
#Peut être utilisée avec csv ou bed
if __name__ == '__main__':
    #Trusted motif list
    firstDF = sys.argv[1]

    #Test motif list
    secondDF = sys.argv[2]
    print(sys.argv[2])
    df1 = pd.read_csv(firstDF, sep='\t')
    df2 = pd.read_csv(secondDF, sep='\t')
    df1 = df1.sort_values("Name")
    df2 = df2.sort_values("Name")
    df1 = df1.reset_index()
    df2 = df2.reset_index()
    #df1 = df1.drop(['coord_to', 'index'], axis = 1)
    #df2 = df2.drop(['coord_to','index'], axis = 1)
    #Imprime la liste des éléments identique entre les deux listes
    df_merge =df1.merge(df2, how="outer",indicator=True, on=['Name','Start',"End"])

    df_merge = df_merge.sort_values(["Name","Start"])
    #df_merge = df_merge.loc[lambda x: x['_merge'] == 'both']
    df_merge.to_csv('DIPPA_BOTH.bed', sep='\t')


    motifChangeList = []
    leftonlyList = []
    rightonlyList = []
    i = 0
    while i < (df_merge.shape[0] - 1):
        row = df_merge.iloc[i,:]
        nextrow= df_merge.iloc[i+1,:]
        name = row["Name"]
        start = row["Start"]
        end = row["End"]
        status = row["_merge"]
        name2 = nextrow["Name"]
        start2 = nextrow["Start"]
        end2 = nextrow["End"]
        status2 = nextrow["_merge"]
        if (status != "both"):
            if status2 != "both" and (name == name2):
                #if (row["_merge"]) == "right_only":
                lengthOverlaps = max(0, min(end,end2) - max(start,start2))
                #interval1 = pd.Interval(left=start, right=end, closed="neither")
                #interval2 = pd.Interval(left=start2, right=end2, closed="neither")

                if (lengthOverlaps > 9):
                    motifChangeList.append([row[["Name","Start","End", "Motif_x", "Motif_y"]].to_list(),nextrow[["Name","Start","End","Motif_x","Motif_y"]].to_list()])
                    i += 1
                else:
                    if status == "right_only":
                        rightonlyList.append(row)
                    if status == "left_only":
                        leftonlyList.append(row)
            else:
                if status == "right_only":
                    rightonlyList.append(row)
                if status == "left_only":
                    leftonlyList.append(row)

        i += 1
    leftmessage = "Left_only Hits: " + str(len(leftonlyList))# + str(leftonlyList)
    rightmessage = "Right_only Hits: " + str(len(rightonlyList))# + str(rightonlyList)
    changemessage = "Motifs changed: " + str(len(motifChangeList))# + str(motifChangeList)
    fulldict = {leftmessage:leftonlyList,rightmessage:rightonlyList,changemessage:motifChangeList}
    print(str(sys.argv[2][:-9]))
    with open(str(sys.argv[2])[:-9] +"/Stats.txt", "w") as fp:
        fp.write(leftmessage + "\n" + rightmessage + "\n" + changemessage + "\n" + str(leftonlyList) + "\n" + str(rightonlyList) + "\n" +  str(motifChangeList))

    #Trouver les éléments qui ne sont pas communs aux deux dataframes
    result1=df_merge.loc[lambda x: x['_merge'] == 'left_only']
    result1=df_merge.loc[lambda x: x['_merge'] == 'left_only']
    #result1 = df1.merge(df2, indicator=True, how='left')#.loc[lambda x: x['_merge'] != 'both']
    result2 = df_merge.loc[lambda x: x['_merge'] == 'right_only']
    #result2 = df1.merge(df2, indicator=True, how='right')#.loc[lambda x: x['_merge'] != 'both']
    #Contient les éléments qui sont seulement dans la première table
    result1.to_csv(str(sys.argv[2][:-9]) +'/DIPPA_COMPARE_LEFT.csv', sep=',')
    result1.to_csv(str(sys.argv[2][:-9]) +'/DIPPA_COMPARE_LEFT.bed', sep='\t')
    # Contient les éléments qui sont seulement dans la deuxième table
    result2.to_csv(str(sys.argv[2][:-9]) +'/DIPPA_COMPARE_RIGHT.csv', sep=',')
    result2.to_csv(str(sys.argv[2][:-9]) +'/DIPPA_COMPARE_RIGHT.bed', sep='\t')