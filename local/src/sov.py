#!/usr/local/bin/python3

'''
Residue-level indexes computed from confusion matrices are not sufficient to properly score a SS prediction method.
Assessing the biological significance of predictions requires comparing predicted and observed segments of residues in the same conformation
• α-Helix segments shorter than 4 residues are meaningless in biology;
• Strands usually involves 2 or more residues (except for isolated β-bridges)
Specific segment-level indexes are then required to analyse these aspects and assess the biological significance of SS predictions.
'''

import pandas as pd
import numpy as np
import argparse

def sov(expected_file, predicted_file):
    secondary_structure = ['H', 'E', '-']
    ss_dictionary = {'H':[0,0,0], 'E':[0,0,0], '-':[0,0,0]}
    with open(expected_file) as filein1, open(predicted_file) as filein2:
        for dssp,pred in zip(filein1, filein2):
            dssp = dssp.rstrip()
            pred = pred.rstrip()
            try: (len(dssp) + len(pred)) % 2 == 0
            except:
                print("Different lenght between expected and predicted sequence: \n%s\n%s" %(dssp,pred))
                raise SystemExit
            else:
                for ss in secondary_structure:
                    num_of_ss = int(dssp.count(ss))
                    dssp_fragments = sov_parser(dssp, ss)
                    pred_fragments = sov_parser(pred, ss)
                    if num_of_ss != 0: 
                        ss_dictionary[ss][0] += sov_scorer(dssp_fragments, pred_fragments, num_of_ss)
                        ss_dictionary[ss][1] += 1
    for val in ss_dictionary.values():
        val[2] = val[0] / val[1]

    return(ss_dictionary)


def sov_parser(sequence, ss):
    '''
    This function is intended to compute the Segment OVerlap index.
    A segment-level index which evaluated the percent average overlap 
    between predicted and observed segments in the three SS conformations.
    '''
    sequence = sequence.rstrip()
    fragments = []
    val = 0
    while val < len(sequence):
        tmp_list = []
        while sequence[val] == ss and val < len(sequence) - 1:
            tmp_list.append(val)
            val += 1
        fragments.append(set(tmp_list))
        val += 1
    
    return fragments


def sov_scorer(dssp_fragments, pred_fragments, num_of_ss):
    summatory = []
    for s in dssp_fragments:
        for t in pred_fragments:
            if s & t:
                minov = len(s&t)
                maxov = len(s|t)
                delta = min([minov, maxov-minov, len(s)/2, len(t)/2])
                summatory.append(((minov + delta)/maxov * len(s)))
    if num_of_ss:
        sov = sum(summatory) * 100 * (1/num_of_ss)
    
    return sov


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