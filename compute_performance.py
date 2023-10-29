
import os
from os import listdir
from os.path import isfile, join
import pandas as pd
import numpy as np
import sys
import sklearn
print(sklearn.__version__)
#from sklearn.metrics import precision_recall_curve
#from sklearn.metrics import roc_curve
#from sklearn.metrics import PrecisionRecallDisplay
#from sklearn.metrics import RocCurveDisplay
#from sklearn.metrics import f1_score
#from sklearn.metrics import auc
from matplotlib import pyplot
from sklearn import metrics


#Permet de comparer deux tables de listes
#python3 compute_performance.py C:\Users\felix\VraiDesktop\Stage-Copie\sortedTrustedPPR_April2023newproteome.txt  C:\Users\felix\VraiDesktop\Stage-Copie\InconclusiveNoAnkyrin.bed  C:\Users\felix\VraiDesktop\Stage-Copie\TPRlist.txt  C:\Users\felix\VraiDesktop\Stage-Copie\Ankyrin.txt BOTHNoCandidateStructureFiltered1.bed Cleaned_NegativeMotifs.txt
def curves(dataframe,mode):
    dataframe = dataframe.dropna(subset=["Label"])
    labels = dataframe["Label"].astype(int)
    evalues = dataframe["E_value"].astype(float).apply(lambda x: -x)
    c_evalues = dataframe["c_Evalue"].astype(float).apply(lambda x: -x)
    i_evalues = dataframe["i_Evalue"].astype(float).apply(lambda x: -x)
    scores = dataframe["score2"].astype(float)
    dataframe.to_csv('example.csv')
    _, ax = pyplot.subplots(figsize=(7, 8))
    _, ax2 = pyplot.subplots(figsize=(7, 8))
    if mode == "PR":
        precisionScore, recallScore, thresholdScore = metrics.precision_recall_curve(labels,scores)
        disp = metrics.PrecisionRecallDisplay(precision=precisionScore, recall=recallScore)
        disp.plot(ax=ax,color='blue',name="Score")
        f1Score = np.array((2*(precisionScore*recallScore))/(precisionScore+recallScore))[:-1]

        print("F1Score: ")
        print(np.mean(f1Score))
        aucScore = metrics.auc(recallScore,precisionScore)
        print("Scores: \nAUC: " + str(np.mean(aucScore)))# + "\nThresholds: " + str(thresholdScore))
        #pyplot.plot(recallScore,precisionScore,color='blue',label="Score")

        precisionEvalue, recallEvalue, thresholdEvalue = metrics.precision_recall_curve(labels, evalues)
        disp = metrics.PrecisionRecallDisplay(precision=precisionEvalue, recall=recallEvalue)
        disp.plot(ax=ax,color='green',name="Evalue")
        f1Evalue = np.array((2 * (precisionEvalue * recallEvalue)) / (precisionEvalue + recallEvalue))[:-1]

        print("F1Evalue: ")
        print(np.mean(f1Evalue))
        #precisionEvalue, recallEvalue, thesholdEvalue = roc_curve(labels, evalues)
        aucEvalue = metrics.auc(recallEvalue, precisionEvalue)
        print("Evalue: \nAUC: " + str(np.mean(aucEvalue)))# + "\nThresholds: " + str(thresholdEvalue))
        #pyplot.plot(recallEvalue, precisionEvalue,color='green',label="Evalue")

        precisionC, recallC, thresholdC = metrics.precision_recall_curve(labels, c_evalues)
        #precisionC, recallC, thesholdC = roc_curve(labels, c_evalues)
        disp = metrics.PrecisionRecallDisplay(precision=precisionC, recall=recallC)
        disp.plot(ax=ax,color='red',name="c-Evalue")
        aucC = metrics.auc(recallC, precisionC)
        f1C = np.array((2 * (precisionC * recallC)) / (precisionC + recallC))[:-1]

        #f1C.plot(ax=ax2, color="lightcoral", name="Score")
        print("F1C: ")
        print(np.mean(f1C))
        print("Conditional: \nAUC: " + str(np.mean(aucC)))#+ "\nThresholds: " + str(thresholdC))
        #pyplot.plot(recallC, precisionC,color='red',label="c-Evalue")

        precisionI, recallI, thresholdI = metrics.precision_recall_curve(labels, i_evalues)
        #precisionI, recallI, thesholdI = roc_curve(labels, i_evalues)
        disp = metrics.PrecisionRecallDisplay(precision=precisionI, recall=recallI)
        disp.plot(ax=ax,color='purple',name="i-Evalue")
        aucI = metrics.auc(recallI, precisionI)
        f1I = np.array((2 * (precisionI * recallI)) / (precisionI + recallI))[:-1]

        #f1I.plot(ax=ax2, color="thistle", name="Score")
        print("F1I: ")
        print(np.mean(f1I))

        print("Independent: \nAUC: " + str(np.mean(aucI)))#+ "\nThresholds: " + str(thresholdI))
        #pyplot.plot(recallI, precisionI,color='purple',label="i-Evalue")
        ax.set_title("Precision-Recall")
        pyplot.xlabel("Thresholds")
        pyplot.ylabel("F1")

        pyplot.subplot(221)
        pyplot.plot(thresholdScore,f1Score, color="lightsteelblue",label="Score")
        pyplot.title("Score")
        pyplot.subplot(222)
        pyplot.plot(-thresholdEvalue,f1Evalue, color="chartreuse",label="Evalue")
        #pyplot.xlim(-0.2,0)
        pyplot.xscale("log")
        pyplot.title("Evalue")
        pyplot.subplot(223)
        pyplot.plot(-thresholdC, f1C, color="lightcoral", label="Conditional")
        pyplot.title("Conditional")
        pyplot.subplot(224)
        pyplot.plot(-thresholdI, f1I, color="thistle", label="Independant")
        pyplot.title("Independant")
        pyplot.subplots_adjust(top=0.92, bottom=0.08, left=0.10, right=0.95, hspace=0.25,wspace=0.35)
        pyplot.show()

    elif mode == "ROC":
        fprScore, tprScore, thesholdScore = metrics.roc_curve(labels, scores)
        aucScore = metrics.auc(fprScore,tprScore)
        print("Score " + str(aucScore))
        display = metrics.RocCurveDisplay(fpr=fprScore,tpr=tprScore,roc_auc=aucScore,estimator_name="Score ROC")
        display.plot(ax=ax,color='blue', label="Score")

        fprE, tprE, thesholdE = metrics.roc_curve(labels, evalues)
        aucE = metrics.auc(fprE,tprE)
        print("Evalue " + str(aucE))
        displayE = metrics.RocCurveDisplay(fpr=fprE,tpr=tprE,roc_auc=aucE,estimator_name="Score ROC")
        displayE.plot(ax=ax,color='green', label="Evalue")

        fprC, tprC, thesholdC = metrics.roc_curve(labels, c_evalues)
        aucC = metrics.auc(fprC,tprC)
        print("cEvalue " + str(aucC))
        displayC = metrics.RocCurveDisplay(fpr=fprC,tpr=tprC,roc_auc=aucC,estimator_name="Score ROC")
        displayC.plot(ax=ax,color='red', label="c-Evalue")

        fprI, tprI, thesholdI = metrics.roc_curve(labels, i_evalues)
        aucI = metrics.auc(fprI,tprI)
        print("iEvalue " + str(aucI))
        displayI = metrics.RocCurveDisplay(fpr=fprI,tpr=tprI,roc_auc=aucI,estimator_name="Score ROC")
        displayI.plot(ax=ax,color='purple', label="i-Evalue")
        ax.set_title("ROC")



def compareTablesV3(df1, df2,label):
    df1["Start"] = df1["Start"].astype(int)
    df2["Start"] = df2["Start"].astype(int)
    df1 = df1.sort_values(by=["Name","Start"])
    df2 = df2.sort_values(by=["Name","Start"])
    df1 = df1.reset_index()
    df2 = df2.reset_index()

    #Imprime la liste des éléments identique entre les deux listes
#    df_merge =df1.merge(df2, how="outer",indicator=True, )

    df_merge =df1.merge(df2, how="outer",indicator=True )
    df_merge = df_merge.sort_values(["Name","Start"])
    #df_merge = df_merge.loc[lambda x: x['_merge'] == 'both']
    df_merge.to_csv('DIPPA_BOTH.bed', sep='\t')


    motifChangeList = []
    leftonlyList = []
    rightonlyList = []
    #For the performance curve, a label of 1 signals a true positive, 0 is a "False positive"
    labelList = []
    i = 0
    newDataframe = pd.DataFrame()
    print(df_merge)
    while i < df_merge.shape[0]:
        row = df_merge.iloc[i,:]
        name = row["Name"]
        start = row["Start"]
        end = row["End"]
        status = row["_merge"]
        eV = row["E_value"]
        cE = row["c_Evalue"]
        iE = row["i_Evalue"]
        newrow = row.copy()
        if i < (df_merge.shape[0]-1):
            nextrow = df_merge.iloc[i + 1, :]
            name2 = nextrow.loc["Name"]
            start2 = nextrow.loc["Start"]
            end2 = nextrow.loc["End"]
            status2 = nextrow.loc["_merge"]


            if (status != "both"):
                if status2 != "both" and (name == name2):
                    #if (row["_merge"]) == "right_only":
                    lengthOverlaps = max(0, min(end,end2) - max(start,start2))

                    if (lengthOverlaps > 9):
                        motifChangeList.append([row[["Name","Start","End", "_merge"]].to_list(),nextrow[["Name","Start","End","_merge"]].to_list()])
                        if status == "right_only":
                            newrow["Label"] = label

                        labelList.append(label)

                        labelList.append(label)
                        i += 1
                        row = df_merge.iloc[i, :]
                        newrow = row.copy()
                        if status2 == "right_only" or status2 =="both":
                            newrow["Label"] = label
                    else:
                        if status == "right_only":
                            rightonlyList.append(row)
                            #newrow["Label"] = np.nan
                        if status == "left_only":
                            leftonlyList.append(row)
                            #newrow["Label"] = np.nan
                else:
                    if status == "right_only":
                        rightonlyList.append(row)
                        #newrow["Label"] = np.nan
                    if status == "left_only":
                        leftonlyList.append(row)
            else:
                #newrow["Label"] = np.nan
                labelList.append(label)
        else:
            # Handles last row
            if status == "right_only":
                rightonlyList.append(row)
                #newrow["Label"] = np.nan
                labelList.append(0)
            elif status == "left_only":
                leftonlyList.append(row)
            elif status == "both":
                newrow["Label"] = label
        if newrow.array[-1] != "left_only":
            newDataframe = pd.concat([newDataframe,newrow],ignore_index=True,axis=1)

        i += 1

    newDataframe = newDataframe.T
    newDataframe = newDataframe[["Name","Start","End","Motif","E_value","c_Evalue","i_Evalue","score2","Label"]]
    #newDataframe.rename(columns={"Motif_y":"Motif"})
    leftmessage = "Left_only Hits: " + str(len(leftonlyList))# + str(leftonlyList)
    rightmessage = "Right_only Hits: " + str(len(rightonlyList))# + str(rightonlyList)
    changemessage = "Motifs changed: " + str(len(motifChangeList))# + str(motifChangeList)
    #fulldict = {leftmessage:leftonlyList,rightmessage:rightonlyList,changemessage:motifChangeList}
    print("Hits missing from trusted list: " + str(len(leftonlyList)))
    print("Total hits from trusted motif list: " + str(df1.shape[0] -len(leftonlyList)))
    print("See specific hits in Stats.txt")
    with open('Performance_Results'+ str(sys.argv[7]) +"/Stats.txt", "w") as fp:
        fp.write(leftmessage + "\n" + rightmessage + "\n" + changemessage + "\n" + str(leftonlyList) + "\n" + str(rightonlyList) + "\n" +  str(motifChangeList))


    #Trouver les éléments qui ne sont pas communs aux deux dataframes
    result1=df_merge.loc[lambda x: x['_merge'] == 'left_only']
    result2 = df_merge.loc[lambda x: x['_merge'] == 'right_only']

    #Contient les éléments qui sont seulement dans la première table
    result1.to_csv('Performance_Results'+ str(sys.argv[7]) +'/DIPPA_COMPARE_LEFT.csv', sep=',')
    result1.to_csv('Performance_Results'+ str(sys.argv[7]) +'/DIPPA_COMPARE_LEFT.bed', sep='\t')
    # Contient les éléments qui sont seulement dans la deuxième table
    result2.to_csv('Performance_Results'+ str(sys.argv[7]) +'/DIPPA_COMPARE_RIGHT.csv', sep=',')
    result2.to_csv('Performance_Results'+ str(sys.argv[7]) +'/DIPPA_COMPARE_RIGHT.bed', sep='\t')


    return newDataframe

if __name__ == '__main__':
    # First DF is the trusted list
    firstDF = sys.argv[1]
    # Second List is the Inconclusive list
    secondDF = sys.argv[2]

    # Third List is the TPR list
    thirdDF = sys.argv[3]

    # Next is the ankyrin list
    ankyrinFile = sys.argv[4]

    #Then, the list of trusted motifs
    motifsFile = sys.argv[5]

    #Then, the list of NEGATIVE motifs
    negativeFile = sys.argv[6]

    # Finally is the list of hits
    hitDF = sys.argv[7]

    directoryName = 'Performance_Results' + str(sys.argv[7])
    os.makedirs(directoryName, exist_ok=True)


    df1 = pd.read_csv(firstDF, sep='\t',index_col=False)
    dfInconclusive = pd.read_csv(secondDF, sep='\t', index_col=False,names=["Name"])
    dfInconclusive = pd.DataFrame(dfInconclusive["Name"])
    dfTPR = pd.read_csv(thirdDF, sep='\t', index_col=False, names=["Name"])
    dfTPR = pd.DataFrame(dfTPR["Name"])
    dfAnkyrin = pd.read_csv(ankyrinFile, sep='\t', index_col=False, names=["Name"])
    dfAnkyrin = pd.DataFrame(dfAnkyrin["Name"])
    dfNegative = pd.read_csv(negativeFile, sep='\t', index_col=False, names=["Name","Start","End"])
    dfPositive = pd.read_csv(motifsFile, sep='\t', index_col=False, header=0)

    #df2 = pd.read_csv(thirdDF, sep='\t',index_col=False,names=["Name", "Start", "End", "Motif"],header=0)
    df2 = pd.read_csv(hitDF, sep='\t', index_col=False, header=0)
    df1 = df1.reset_index()
    dfInconclusive = dfInconclusive.reset_index()
    dfTPR = dfTPR.reset_index()
    dfAnkyrin = dfAnkyrin.reset_index()
    df2 = df2.reset_index()
    #Imprime la liste des éléments identique entre les listes
    df_mergeTrust =df1.merge(df2, how="inner",on="Name")
    df_mergeTrust=df_mergeTrust.drop(axis=1, labels=["index_x","index_y"])

    df_mergeTrust.to_csv('Performance_Results'+ str(sys.argv[7]) + '/BOTH.csv', sep=',')
    df_mergeTrust.to_csv('Performance_Results'+ str(sys.argv[7]) + '/BOTH.bed', sep='\t',index=False)


    #Trouver les éléments qui ne sont pas communs aux deux dataframes
    df2names = df2[["Name","tlen","Start","End","Motif","E_value","c_Evalue","i_Evalue","score2"]]
    #The sequences only in the trusted file and not in the result file
    result1 = df1.merge(df2names, indicator=True, how='left',on="Name").loc[lambda x: x['_merge'] != 'both']
    #The hits in the result file and not the trusted file
    result2 = df1.merge(df2names, indicator=True, how='right',on="Name").loc[lambda x: x['_merge'] != 'both']
    result2Names = pd.DataFrame(result2[["Name","tlen","Start","End","Motif","score2"]])

    #The hits in the resultfile and not in the trusted OR inconclusive file (Should still contain TPR and Ankryin + unknowns)
    result3NoInc = result2Names.merge(dfInconclusive,indicator=True,how='left',on="Name").loc[lambda x: x['_merge'] != 'both']
    result3NoInc = result3NoInc.drop(["_merge"], axis=1)
    # The hits in the resultfile and not in the trusted NOR TPR file NOR Inconclusive (Should still contain Ankyrin + unknowns)
    result3NoIncTPR = result3NoInc.merge(dfTPR, indicator=True, how='left',on="Name").loc[lambda x: x['_merge'] != 'both']
    result3NoIncTPR = result3NoIncTPR.drop(["_merge"], axis=1)

    #The hits in the resultfile and not in the trusted NOR TPR file NOR Inconclusive, NOR Ankyrin (rest is unknown)
    result3Unknowns = result3NoIncTPR.merge(dfAnkyrin, indicator=True, how='left',on="Name").loc[lambda x: x['_merge'] != 'both']
    result3Unknowns = result3Unknowns.drop(["_merge"], axis=1)
    result3Ank = result3NoIncTPR.merge(result3Unknowns, indicator=True, how='left',on="Name").loc[lambda x: x['_merge'] != 'both']

    # Contains only TPR + unkowns
    result3NoIncAnk = result3NoInc.merge(dfAnkyrin, indicator=True, how='left',on="Name").loc[lambda x: x['_merge'] != 'both']
    result3NoIncAnk = result3NoIncAnk.drop(["_merge"], axis=1)
    result3NoIncAnk = result3NoIncAnk.merge(result3Unknowns, indicator=True, how='left',on="Name").loc[lambda x: x['_merge'] != 'both']

    # Contains no TPR
    result2.drop(["_merge"],axis=1,inplace=True)
    result3NoTPR = result2.merge(dfTPR, indicator=True, how='left',on="Name").loc[lambda x: x['_merge'] != 'both']
    result3NoTPR = result3NoTPR.drop(["_merge"], axis=1)

    # Contains no TPR and Ank (Inc + unknowns)
    result3NoTPRAnk = result3NoTPR.merge(dfAnkyrin, indicator=True, how='left',on="Name").loc[lambda x: x['_merge'] != 'both']
    result3NoTPRAnk = result3NoTPRAnk.drop(["_merge"], axis=1)
    result3NoTPRAnk = result3NoTPRAnk.merge(result3Unknowns, indicator=True, how='left',on="Name").loc[lambda x: x['_merge'] != 'both']


    #result2 = result2.drop(["_merge"],axis=1) #Not trusted

    result3Inc = result3NoTPRAnk.drop(["index_x_x"], axis=1) #Inconclusives only
    result3Inc = result3Inc.iloc[:,0:5]
    result3Inc = result3Inc.rename(columns={"Name":"Name","tlen_x":"tlen","Start_x":"Start","End_x":"End","Motif_x":"Motif","E_value_x":"E_value","c_Evalue_x":"c_Evalue","i_Evalue_x":"i_Evalue","score2_x":"score2"})
    result3TPR = result3NoIncAnk.drop(["_merge"], axis=1) # TPR only
    result3TPR = result3TPR.iloc[:, 0:5]
    result3TPR = result3TPR.rename(columns={"Name": "Name", "tlen_x": "tlen", "Start_x": "Start", "End_x": "End", "Motif_x": "Motif","E_value_x":"E_value","c_Evalue_x":"c_Evalue","i_Evalue_x":"i_Evalue","score2_x":"score2"})
    result3Ank = result3Ank.drop(["_merge"], axis=1) # Ank only
    result3Ank = result3Ank.iloc[:, 0:5]
    result3Ank = result3Ank.rename(columns={"Name": "Name", "tlen_x": "tlen", "Start_x": "Start", "End_x": "End", "Motif_x": "Motif","E_value_x":"E_value","c_Evalue_x":"c_Evalue","i_Evalue_x":"i_Evalue","score2_x":"score2"})


    print("Number of Trusted hits:" + str(df_mergeTrust.shape[0]))
    motif_mergeSeries = df_mergeTrust["Motif"].value_counts(normalize=False)
    print(motif_mergeSeries)


    motif_mergeSeries.to_csv('MotifTrusted.csv', sep=',')
    print("Number of Inconclusive hits:" + str(result3Inc.shape[0]))
    motif_mergeSeries2Inc = result3Inc["Motif"].value_counts(normalize=False)
    print(motif_mergeSeries2Inc)
    print("Number of TPR hits:" + str(result3TPR.shape[0]))
    motif_mergeSeries2TPR = result3TPR["Motif"].value_counts(normalize=False)
    print(motif_mergeSeries2TPR)

    print("Number of Ankyrin hits:" + str(result3Ank.shape[0]))
    motif_mergeSeries2Ank = result3Ank["Motif"].value_counts(normalize=False)
    print(motif_mergeSeries2Ank)

    print("Number of Unknowns:" + str(result3Unknowns.shape[0]))
    motif_mergeSeries3 = result3Unknowns["Motif"].value_counts(normalize=False)
    print(motif_mergeSeries3)

    print("Total hits: " + str(
        df_mergeTrust.shape[0] + result3Inc.shape[0] + result3TPR.shape[0] + result3Ank.shape[0] +
        result3Unknowns.shape[0]))

    #Contient les éléments qui sont seulement dans la première table
    result1.to_csv('Performance_Results'+ str(sys.argv[7]) + '/MissingFromTrusted.bed', sep='\t',index=False)
    # Contient les éléments qui sont seulement dans la deuxième table
    result3Unknowns.to_csv('Performance_Results'+ str(sys.argv[7]) + '/Unknowns.bed', sep='\t')
    result3Inc.to_csv('Performance_Results'+ str(sys.argv[7]) + '/COMPARE_Inconclusive.bed', sep='\t',index=False)
    result3Ank.to_csv('Performance_Results'+ str(sys.argv[7]) + '/COMPARE_Ankyrin.bed', sep='\t', index=False)
    result3TPR.to_csv('Performance_Results'+ str(sys.argv[7]) + '/COMPARE_TPR.bed', sep='\t', index=False)

    df2NoCandidate = df2.loc[df2["Motif"]!="Candidate++"]
    print("Band")
    print(df2NoCandidate)
    df2NoCandidate["Labels"] = np.nan
    trustedlabeledList = compareTablesV3(dfPositive, df2NoCandidate,1)
    print("trustedlabaledlist")
    print(trustedlabeledList)
    fulllabeledList = compareTablesV3(dfNegative, trustedlabeledList, 0)
    #fulllabeledList = fulllabeledList.dropna()
    # Contient la liste originale avec les labels
    fulllabeledList.to_csv('Performance_Results'+ str(sys.argv[7]) +'/LabeledList.csv', sep=',')
    fulllabeledList.to_csv('Performance_Results'+ str(sys.argv[7]) +'/LabeledList.bed', sep='\t')
    curves(fulllabeledList,"PR")