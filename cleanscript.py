import pandas as pd
import sys


if __name__ == '__main__':
    file = sys.argv[1]
    with open("TabCleaned" + file,'w') as fw:
        with open(file, 'r') as fr:
            lines = fr.readlines()
            print(len(lines))
            for line in lines:
                print(line)
                while line[-2] == "\t":
                    line = line[:-2] + "\n"
                    print(line)
                #line.replace("\t\t","\t")
                #line.replace("\t\n","\n")
                fw.write(line)

    dataframe = pd.read_csv("TabCleaned" + file,sep="\t",names=["Name","tlen","Start","End","Motif","Extra"])
    print(dataframe.columns)
    dataframe2 = dataframe[["Name","Start","End"]]
    dataframe2.to_csv("Cleaned_" + file,sep="\t",index=False)
