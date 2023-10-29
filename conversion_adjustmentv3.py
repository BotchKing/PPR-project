#!/usr/bin/env python3
import sys

import numpy as np
import csv
from re import split
import pandas as pd


def hmm_adjust_coord(converted):
    df1 = pd.read_csv(converted + "_converted.csv",sep="\t")
    df2 = df1.copy()
    # Remove sequences if the score is less than 30
    # we don't want junk sequences when we build new profiles, it would reduce specificity
    #df2 = df2.loc[df2["score2"] >= 10]

    #remove targets that length is less than profiles HMM's
    df2 = df2.loc[df2["tlen"] >= df2["qlen"]]
    # adjusted targets length
    # HMMER is 1-based indexing (count starts 1...N for sequence of N length)
    # Sequence provided by proteom is 0-based indexing
    # So these coordinates below are 0-based indexing
    # More details: http://genome.ucsc.edu/blog/the-ucsc-genome-browser-coordinate-counting-systems/
    df2 = df2.loc[df2["tlen"] >= df2["qlen"]]
    df2["ali_from"] = df2["ali_from"] - (df2["hmm_from"]) + 1
    #df2["ali_from"].apply(lambda x: 0 if x < 0)
    df2["ali_to"] = df2["ali_to"] + (df2["qlen"] - df2["hmm_to"])

    df2 = df2[df2['ali_from'] > 0]
    df2 = df2[df2['ali_to'] <= df2['tlen']]
    df2 = df2[["target_name", "tlen", "query_name", "qlen", "E_value",
               "c_Evalue", "i_Evalue", "ali_from", "ali_to","score2"]]
    #
    df2.to_csv((str(filename) + "_adjusted.csv"),sep = '\t', index = None, float_format = '%g')
    print('Targets coordinates adjusted...')


def hmmsearch_output_to_csv(filename1):

        hmm_hits = []
        header = ["target_name", "accession", "tlen", "query_name", "accession", "qlen",
                  "E_value", "score", "bias", "#", "of", "c_Evalue", "i_Evalue", "score2",
                  "bias", "hmm_from", "hmm_to", "ali_from", "ali_to", "env_from", "env_to", 
                  "acc", "description_of_target"]

        hmm_hits.append(header)
    
        with open (filename1,'r') as file :
            f_list = [line.rstrip('\n') for line in file]
            for line in f_list:
                if line.startswith("#") or line == "":
                     continue
                line = split("\s+",line, 22)
                hmm_hits.append(line)
        hits_csv = np.asarray(hmm_hits)
        np.savetxt(str(filename1)+"_converted.csv", hits_csv, delimiter = "\t", fmt = '%s')
        print('conversion is done...')
        return (hits_csv)
        
if __name__ == '__main__':
    filename = sys.argv[1]
    converted_file = hmmsearch_output_to_csv(filename)
    hmm_adjust_coord(filename)
           

