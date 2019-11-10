#!/usr/local/bin/python3
# import sys
from sklearn.metrics import multilabel_confusion_matrix
import pandas as pd
import numpy as np
import argparse

def pretty_dictionary(D):
    for key,val in D.items():
        print(key, " = ", val)
    print()

def compute_similarity(expected, predicted):
    # with open(expected) as dssp, open(predicted) as predic:
    true = []
    pred = []
    for i,j in zip(expected, predicted):
        i = i.rstrip()
        j = j.rstrip()
        for k in range(len(i)):
            true.append(i[k])
            pred.append(j[k])
    return(true, pred)

#####################################
### One-vs-all confusion matrices ###
#####################################

def print_performance(confusion_matrix):
    TP = confusion_matrix[1][1]
    FP = confusion_matrix[0][1]
    TN = confusion_matrix[0][0]
    FN = confusion_matrix[1][0]
        
    ACC = (TP + TN)/(TP + FP + TN + FN)
    TPR = TP / (TP + FN)
    TNR = TN / (TN + FP)
    FPR = 1 - TNR
    FNR = 1 - TPR
    PPV = TP / (TP + FP)
    NPV = TN / (TN + FN)
                    
    MCC = ((TP * TN) - (FN * FP)) / (np.sqrt((TP + FN) * (TP + FP) * (TN + FN) * (TN + FP)))
   
    return(MCC, ACC, TPR, PPV, FPR, NPV)


def ss_compo(dssp_file):
    ss = ['H', 'E', '-']
    dictionary={'H':0, 'E':0, '-':0}
    for line in dssp_file:
        line = str(line.rstrip())
        for ch in line:
            dictionary[ch] += 1

    return (dictionary)
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('DSSP', \
        help='This is the dssp file with concatenating BlindSet secondary structure')
    parser.add_argument('SSpredictions', \
        help='This is the file with concatenating BlindSet secondary structure predicted via GOR')
    args = parser.parse_args()    
    
    #ss_composition = ss_compo(dssp)
    #pretty_dictionary(ss_composition)
    
    x,y = compute_similarity(args.DSSP, args.SSpredictions)
    multiclass_conf_mat = multilabel_confusion_matrix(x,y, labels=['H','E','-'])
    total = np.sum(multiclass_conf_mat[0])
    print('\nHelix one-vs-all:\n', multiclass_conf_mat[0], '\n\nStrand one-vs-all:\n', multiclass_conf_mat[1], '\n\nCoil one-vs-all:\n', multiclass_conf_mat[2] )
    ss = ['H', 'E', 'C']
    true_predicted = 0

    dictionary = {'H': None, 'E': None, 'C': None}
    for index in range(3):
        true_predicted += multiclass_conf_mat[index][1][1]
        MCC, ACC, TPR, PPV, FPR, NPV = print_performance(multiclass_conf_mat[index])
        dictionary[ss[index]] = [MCC, ACC, TPR, PPV, FPR, NPV]
    d = pd.DataFrame(data=dictionary, index=['MCC', 'ACC', 'TPR', 'PPV', 'FPR', 'NPV'], dtype=float)

    q3_score = true_predicted/total
    print(d)
    print('\nThis is the Q3 score:\t%f\n' %q3_score)