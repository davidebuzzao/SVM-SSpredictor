#!/usr/local/bin/python3
import sys, pickle, pandas as pd, numpy as np
from Svm import Svm
from Tools import Dataset, Pssm, Dssp
from sov import sov, sov_parser, sov_scorer
from performance import compute_similarity, print_performance
from sklearn.metrics import multilabel_confusion_matrix

'''
    Input file:
    
    Fold               Dataset.pkl                Dataset IDs                Dataset.dat                Dataset.pred
    fold1       dataset/cv/fold1/cv5.pkl    dataset/cv/fold1/cv5.id     dataset/cv/fold1/cv5.dat    dataset/cv/fold1/cv5.pred
'''

if __name__ == '__main__':
    filein = sys.argv[1]
    
    ss = ['H', 'E', '-']
    cv_list = ['fold1', 'fold2', 'fold3', 'fold4', 'fold5']
    scores = ['MCC', 'ACC', 'TPR', 'PPV', 'FPR', 'NPV', 'SOV']
    cv_dictionary = dict([key1, dict([key2, np.zeros(len(ss))] for key2 in scores)] for key1 in cv_list)
    q3_scores = np.zeros(len(cv_list))

    with open(filein) as f:

        for line in f:
            items = line.split()
            print(items)
            
            fold = items[0]
            data_pkl = items[1]
            data_id = items[2]
            data_dat = items[3]
            pred_file = items[4]
            
            ### Build the dataset from dataset.pkl ###
            model = Svm() \
                        .load(dataset=data_pkl, id_file=data_id, encode=True, pkl=True)\
                        .fetch_prediction(prediction_file=pred_file)\
                        #.save(path=data_dat, format='dat')

            dictionary = model.fetch_dictionary()

            '''
            Performance Steps:
            1. Confusion Matrix and related scores
            2. Segment OVerlapping score
            '''
            prediction = []
            expectation = []

            for key in dictionary:
                expectation.append(dictionary[key]['dssp'])
                prediction.append(dictionary[key]['svm_pred'])

            x,y = compute_similarity(expectation, prediction)
            multiclass_conf_mat = multilabel_confusion_matrix(x,y, labels=['H','E','-'])
            total = np.sum(multiclass_conf_mat[0])

            true_predicted = 0

            for index in range(3):
                true_predicted += multiclass_conf_mat[index][1][1]
                MCC, ACC, TPR, PPV, FPR, NPV = print_performance(multiclass_conf_mat[index])
                score = [MCC, ACC, TPR, PPV, FPR, NPV]
                
                scr = 0
                while scr < len(score):
                    cv_dictionary[fold][scores[scr]][index] = round(score[scr]*100,2)
                    scr += 1

            q3 = true_predicted/total
            q3_scores[cv_list.index(fold)] = q3

            sov_scores = sov(expectation, prediction, file=False)
            final_sov = np.zeros(len(ss))
            for i in range(len(ss)):
                final_sov[i] = sov_scores[ss[i]][2]
            
            cv_dictionary[fold]['SOV'] += np.round(final_sov,2)
    
    
    ### Find the way to print better c !!! ###
    d = pd.DataFrame(cv_dictionary)
    
    print('\n',d)
    cv_average_ss = dict([key, np.zeros(len(ss))] for key in scores)
    cv_average = dict([key, np.zeros(1)] for key in scores)

    divisor1 = len(cv_list)
    
    for i in range(len(ss)):
        scr = 0
        while scr < len(scores):
            for j in cv_list:
                cv_average_ss[scores[scr]][i] += cv_dictionary[j][scores[scr]][i]
            cv_average_ss[scores[scr]][i] #/= divisor1
            scr += 1
    d1 = pd.DataFrame(cv_average_ss, index=['H','E','C'])
    print('\n', d1, '\n')
    
    divisor2 = len(ss)
    for i in range(len(scores)):
        cv_average[scores[i]] += np.sum(cv_average_ss[scores[i]])
        cv_average[scores[i]] /= divisor2
    d2 = pd.DataFrame(cv_average, index=['Average values'])
    print(d2, '\n')