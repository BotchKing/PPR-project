#!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np


def cutoff(filename,threshold):

    df1 = pd.read_csv(filename, sep = '\t')
    df2 = df1.copy()
    # Remove sequences if the score is less than 30
    # we don't want junk sequences when we build new profiles, it would reduce specificity

    #removal of sequences that belong in a low Evalue sequence
    df3 = df2.copy()
    grouped = df3["target_name"].unique()

    for names in grouped:
        group = df2[df2["target_name"]==names]
        name = group["target_name"].iloc[0]
        Evalue= group["E_value"]
        minimumEvalue = Evalue.min()
        if minimumEvalue > threshold:
            df3 = df3[df3["target_name"] != name]
            print("Removed " + name + " ; with lowest Evalue at " + str(minimumEvalue))
    df2 = df3

    df2.to_csv(str(filename)[:-4] + "_Evalue_" + str(threshold) + ".csv",sep = '\t', index = None, float_format = '%g')
    print('Removed motifs with lowest sequence Evalue > ' + str(threshold))
if __name__ == '__main__':
    filename = sys.argv[1]

    if len(sys.argv) == 3:
        threshold = int(sys.argv[2])

    else:
        threshold = 0.05

    cutoff(filename,threshold)