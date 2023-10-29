#!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np


def cutoff(type,filename,threshold):
    df1 = pd.read_csv(filename, sep = '\t')
    df2 = df1.copy()

    #score cutoff
    if type == "s":
        if isinstance(threshold, str):
            #default threshold
            threshold = 0
            print("Default score threshold of " + str(threshold) + " used")
        df3 = df2.loc[df2["score2"] >= threshold]
        dfOut = df2.loc[df2["score2"] < threshold]
        print(dfOut)
        print("Removed " + str(dfOut.shape[0]) + " sequences with score under " + str(threshold))
        print("Removed motifs can be found in 'Removed.csv'")
        dfOut.to_csv("Removed.csv", sep=',', index=None, float_format='%g')
        df3.to_csv(str(filename)[:-4] + "_Score_" + str(threshold) + ".csv",sep = '\t', index = None, float_format = '%g')

    if type == "c":
        if isinstance(threshold, str):
            #default threshold
            threshold = 0.05
            print("Default evalue threshold of " + str(threshold) + " used")
        df3 = df2.loc[df2["c_Evalue"] <= threshold]
        dfOut = df2.loc[df2["c_Evalue"] > threshold]
        print(dfOut)
        print("Removed " + str(dfOut.shape[0]) + " motifs with c-Evalue over " + str(threshold))
        print("Removed motifs can be found in 'Removed.csv'")
        dfOut.to_csv("Removed.csv", sep=',', index=None, float_format='%g')
        df3.to_csv(str(filename)[:-4] + "_cEvalue_" + str(threshold) + ".csv", sep='\t', index=None, float_format='%g')

    if type == "i":
        if isinstance(threshold, str):
            #default threshold
            threshold = 0.05
            print("Default evalue threshold of " + str(threshold) + " used")
        df3 = df2.loc[df2["i_Evalue"] <= threshold]
        dfOut = df2.loc[df2["i_Evalue"] > threshold]
        print(dfOut)
        print("Removed " + str(dfOut.shape[0]) + " motifs with i-Evalue over " + str(threshold))
        print("Removed motifs can be found in 'Removed.csv'")
        dfOut.to_csv("Removed.csv", sep=',', index=None, float_format='%g')
        df3.to_csv(str(filename)[:-4] + "_iEvalue_" + str(threshold) + ".csv", sep='\t', index=None, float_format='%g')

    if type == "seq":
        #removed_sequences = np.array([],dtype=object)
        removed_sequences = pd.DataFrame()
        if isinstance(threshold, str):
            #default threshold
            threshold = 0.05
            print("Default evalue threshold of " + str(threshold) + " used")

        grouped = df2["target_name"].unique()
        for names in grouped:
            group = df2[df2["target_name"] == names]
            name = group["target_name"].iloc[0]
            Evalue = group["E_value"]
            minimumEvalue = Evalue.min()
            if minimumEvalue > threshold:
                removed_sequences = pd.concat([removed_sequences,df2.loc[df2["target_name"] == name]])
                df3 = df2.loc[df2["target_name"] != name]
                print("Removed " + name + " ; with lowest Evalue at " + str(minimumEvalue))
        dfOut = pd.DataFrame(removed_sequences)
        print(dfOut)
        print("Removed " + str(dfOut.shape[0]) + " motifs with sequence Evalue over " + str(threshold))
        print("Removed motifs can be found in 'Removed.csv'")
        dfOut.to_csv("Removed.csv", sep=',', index=None, float_format='%g')
        df3.to_csv(str(filename)[:-4] + "_Evalue_" + str(threshold) + ".csv", sep='\t', index=None, float_format='%g')


if __name__ == '__main__':
    filename = sys.argv[1]
    if len(sys.argv) == 2:
        print("Specify the cutoff type. Either 's' for score, 'c' for conditional e-value, 'i' for independant e-value or 'seq' for sequence e-value")
        print("You can specify threshold value, if no value is given a default one will be used (0 for score and 0.05 for e-value-based thresholds)")

    else:
        type = sys.argv[2]
        if len(sys.argv) == 4:
            threshold = int(sys.argv[3])

        else:
            threshold = "default"

        cutoff(type,filename,threshold)