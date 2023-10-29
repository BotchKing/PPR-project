#!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np
import csv
import os

os.makedirs('BED_files', exist_ok= True)

def expandCoords(dataframe1):
    dataframe = dataframe1.copy()
    start2list = []
    end2list = []
    #print(dataframe[dataframe["Motif"] == "E1"])
    #dataframe = dataframe.tail(-1)
    for i in range(dataframe.shape[0]):
        row = dataframe.iloc[i,:]
        if row["Start"] != "Start":
            start = int(row["Start"])
            end = int(row["End"])
            tlen = int(row["tlen"])
            if row["Motif"] == "E1":
                print(row)
                #print(True)
                print(start)

            if start >9:
                start2 = start - 10
            else:
                start2 = 0
            if tlen > (end +10):
                end2 = end + 10
            else:
                end2 = tlen-1
            newrow = row.copy()
            newrow["Start"] = start2
            newrow["End"] = end2

            dataframe.iloc[i,:] = newrow
        #dataframe.iloc[i,:] = row
        #start2list.append(start2)
        #end2list.append(end2)

    #dataframe["Start"] = start2list
    #dataframe["End"] = end2list
    #print(dataframe[dataframe["Motif"] == "E1"])
    return dataframe




#Separates the csv into multiple bed files for each motif, in a new folder
def profile_in_bed():
     filename = sys.argv[1]
     df = pd.read_csv(filename, sep='\t',names=["Name", "tlen", "Motif", "qlen", "E_value",
               "c_Evalue", "i_Evalue", "Start", "End","score2"])

     profile_name = df['Motif'].unique()

     for i in profile_name :
        if i != "Motif":
            dfBED = df[df['Motif'] == i]
            print(i)
            #dfBED.to_csv('./BED_files/'+ i +'Tlen+Score.bed', sep='\t', header = False, index = False)
            dfBEDCut = dfBED[['Name','Start','End']]
            dfBEDCut.to_csv('./BED_files/'+ i +'.bed', sep='\t', header = False, index = False)

            dfBEDexpanded = expandCoords(dfBED)
            #dfBEDexpanded.to_csv('./BED_files/'+ i +'ExpandedTlen+Score.bed', sep='\t', header = False, index = False)
            dfBEDCutExpanded = dfBEDexpanded[['Name','Start','End']]
            dfBEDCutExpanded.to_csv('./BED_files/'+ i +'Expanded.bed', sep='\t', header = False, index = False)
     print('Profiles coordinates are available in this folder: ./BED_files')

if __name__ == '__main__' :
    profile_in_bed()
