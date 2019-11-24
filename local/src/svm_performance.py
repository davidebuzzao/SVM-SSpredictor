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
    try:
        filein = sys.argv[1]
        perftype = sys.argv[2]
    except:
        print('Program Usage: svm_performance <cv_pred.txt> ')
        raise SystemExit
    else:
        num_fold = 5
        ss = ['-', 'H', 'E']
        if perftype == 'cv':
            set_list = ['fold'+ str(i) for i in range(1, num_fold + 1)]
        else:
            set_list = [perftype]

        scores = ['MCC', 'ACC', 'TPR', 'PPV', 'SOV']
        cv_dictionary = dict([key1, dict([key2, np.zeros(len(ss))] for key2 in scores)] for key1 in set_list)
        q3_scores = np.zeros(len(set_list))

        print('\nPerformance is being computed:')
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

                ##############################################
                ## Performance Steps:                       ##
                ## 1. Confusion Matrix and related scores   ##
                ## 2. Segment OVerlapping score             ##
                ##############################################

                prediction = []
                expectation = []
                
                for key in dictionary:
                    expectation.append(dictionary[key]['dssp'])
                    prediction.append(dictionary[key]['svm_pred'])

                x,y = compute_similarity(expectation, prediction)
                multiclass_conf_mat = multilabel_confusion_matrix(x,y, labels=['-','H','E'])
                total = np.sum(multiclass_conf_mat[0])

                true_predicted = 0

                for index in range(len(ss)):
                    true_predicted += multiclass_conf_mat[index][1][1]
                    MCC, ACC, TPR, PPV = print_performance(multiclass_conf_mat[index])
                    score = [MCC, ACC, TPR, PPV]
                    
                    scr = 0
                    while scr < len(score):
                        cv_dictionary[fold][scores[scr]][index] = round(score[scr]*100,2)
                        scr += 1

                q3 = true_predicted/total
                q3_scores[set_list.index(fold)] = round(q3 * 100, 2)

                sov_scores = sov(expectation, prediction)
                final_sov = np.zeros(len(ss))
                for i in range(len(ss)):
                    final_sov[i] = sov_scores[ss[i]][2]
                
                cv_dictionary[fold]['SOV'] += np.round(final_sov,2)

        print('\nPerformance is being showed:')
        df_dict = pd.DataFrame(cv_dictionary)
        if perftype == 'cv':
            print('\nFolds DataFrame:\n\n', df_dict)

        cv_average_ss = dict([key, np.zeros(len(ss))] for key in scores)
        cv_average = dict([key, np.zeros(1)] for key in scores)
        divisor1 = len(set_list)
        for i in range(len(ss)):
            scr = 0
            while scr < len(scores):
                for j in set_list:
                    cv_average_ss[scores[scr]][i] += cv_dictionary[j][scores[scr]][i]
                cv_average_ss[scores[scr]][i] /= divisor1
                scr += 1

        df_meanSS = pd.DataFrame(cv_average_ss, index=['C','H','E'])
        print('\nFolds averaged DataFrame per SS:\n\n', df_meanSS, '\n')

        if perftype == 'cv': 
            cv_stderr_ss = dict([key, np.zeros(len(ss))] for key in scores)
            for i in range(len(ss)):
                scr = 0
                while scr < len(scores):
                    for j in set_list:
                        cv_stderr_ss[scores[scr]][i] += np.power(cv_dictionary[j][scores[scr]][i] - cv_average_ss[scores[scr]][i], 2)
                    
                    cv_stderr_ss[scores[scr]][i] = (np.sqrt(cv_stderr_ss[scores[scr]][i]/(len(set_list)-1)))/np.sqrt(len(set_list))
                    scr += 1
            df_stderr = pd.DataFrame(cv_stderr_ss, index=['C','H','E'])
            print('\nStandard Errors DataFrame:\n\n', df_stderr, '\n')


        divisor2 = len(ss)
        for i in range(len(scores)):
            cv_average[scores[i]] += np.sum(cv_average_ss[scores[i]])
            cv_average[scores[i]] /= divisor2

        df_mean = pd.DataFrame(cv_average, index=['Average values'])
        print('\nFolds averaged DataFrame:\n\n', df_mean, '\n')
        
        q3_mean = np.mean(q3_scores)
        if perftype == 'cv':
            q3_dev = 0
            for i in range(len(q3_scores)):
                q3_dev += np.power(q3_scores[i] - q3_mean, 2)

            q3_stderr = (np.sqrt(q3_dev/(len(set_list)-1)))/np.sqrt(len(set_list))
            print('Q3 scores: ', q3_scores, '\nQ3 mean: ', q3_mean, '±', round(q3_stderr,2))
        else:
            print('Q3 scores: ', q3_scores, '\nQ3 mean: ', q3_mean)