#!/usr/local/bin/python3

import pandas as pd
import numpy as np
import argparse

'''
Residue-level indexes computed from confusion matrices are not sufficient to properly score a SS prediction method.
Assessing the biological significance of predictions requires comparing predicted and observed segments of residues in the same conformation
• α-Helix segments shorter than 3-4 residues are meaningless in biology;
• Strands usually involves 2 or more residues (except for isolated β-bridges)
Specific segment-level indexes are then required to analyse these aspects and assess the biological significance of SS predictions.
'''

def sov(expected, predicted):
    '''
    This function is intended to give back a dictionary of SOV scores, one per secondary structure state. 
    '''
    secondary_structure = ['-', 'H', 'E']
    ss_dictionary = dict([key, np.zeros(3)] for key in secondary_structure)
    
    for dssp,pred in zip(expected, predicted):
            dssp = dssp.rstrip()
            pred = pred.rstrip()
            try: (len(dssp) + len(pred)) % 2 == 0
            except:
                print("Different lenght between expected and predicted sequence: \n%s\n%s" %(dssp,pred))
                raise SystemExit
            else:
                for ss in secondary_structure:
                    dssp_fragments = sov_parser(dssp, ss)
                    pred_fragments = sov_parser(pred, ss)
                    if dssp.count(ss) != 0: 
                        ss_dictionary[ss][0] += sov_scorer(dssp_fragments, pred_fragments)
                        ss_dictionary[ss][1] += 1
    
    for val in ss_dictionary.values():
        val[2] = val[0] / val[1]

    return(ss_dictionary)


def sov_parser(sequence, ss):
    '''
    This function is intended to extract from sequence in input every segment of 
    continous characters of specified secondary structure state
    '''
    sequence = sequence.rstrip()
    fragments = []
    val = 0
    while val < len(sequence):
        tmp_list = []
        while sequence[val] == ss and val < len(sequence) - 1:
            tmp_list.append(val)
            val += 1
        
        if tmp_list and val == len(sequence)-1:     
            tmp_list.append(val)
            
        if len(tmp_list) > 0:
            fragments.append(set(tmp_list))
        val += 1

    return(fragments)


def sov_scorer(dssp_fragments, pred_fragments):
    '''
    This function is intended to compute the Segment OVerlap index (SOV).
    The normalizer factor is computed as introduced by Zemla et al., 1999
    '''
    summatory = []
    normalizer = 0
   
    for obs in dssp_fragments:
        flag = 0
        for pred in pred_fragments:
            if obs & pred:
                flag = 1
                normalizer += len(obs)
                minov = len(obs & pred)
                maxov = len(obs | pred)
                delta = min([minov, maxov-minov, len(obs)//2, len(pred)//2])
                summatory.append(((minov + delta)/maxov * len(obs)))
        if flag == 0:
            normalizer += len(obs)

    if normalizer:
        sov = sum(summatory) * 100 * (1/normalizer)
    else:
        sov = 0

    return(sov)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('DSSP', \
        help='This is the dssp file with concatenating BlindSet secondary structure')
    parser.add_argument('SSpredictions', \
        help='This is the file with concatenating BlindSet secondary structure predicted via GOR')
    args = parser.parse_args()    
    
    sov_scores = sov(args.DSSP, args.SSpredictions)
    dataframe = pd.DataFrame(sov_scores, index=['Sum of SOV values', 'Number of sequences counted', 'SOV value'])
    print(dataframe.to_string(formatters={
                                        'H':'{:,.2f}'.format, 
                                        'E':'{:,.2f}'.format, 
                                        '-':'{:,.2f}'.format}))