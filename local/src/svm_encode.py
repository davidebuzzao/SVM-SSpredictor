#!/usr/local/bin/python3

import sys, pickle, argparse, numpy as np
from Svm import Svm
from Tools import Dataset, Pssm, Dssp

'''
    Input file:
    
    Fold               Dataset.pkl                Dataset IDs                Dataset.dat                Dataset.pred
    fold1       dataset/cv/fold1/cv5.pkl    dataset/cv/fold1/cv5.id     dataset/cv/fold1/cv5.dat    dataset/cv/fold1/cv5.pred
'''

if __name__ == '__main__':
    try:
        filein = sys.argv[1]
        setype = sys.argv[2]
    except:
        print('Program Usage: python3 svm_encode.py <file.txt> <set type (trainingset, testset)>')
        raise SystemExit
    else:
        with open(filein) as f:
            for line in f:
                items = line.split()
                print(items)
                
                path = items[0]
                data_pkl = items[1]
                data_id = items[2]
                data_dat = items[3]

        # Build the dataset from scratch
                prof = Pssm(data_id, setype=setype, raw_file=False).parse()
                dict_prof = prof.fetch_dict()
                dssp = Dssp(data_id, setype=setype, raw_file=False).parse()
                dict_dssp = dssp.fetch_dict()
                dataset = Dataset(data_id, setype=setype).build(profile=dict_prof, dssp=dict_dssp).fetch_dict()

                model = Svm(id_file=data_id, setype=setype)\
                            .load(dataset=dataset, id_file=data_id, encode=True)\
                            .save(path=data_dat, format='dat') \
                            .save(path=data_pkl, format='pkl')