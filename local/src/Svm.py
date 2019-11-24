#!/usr/local/bin/python3
import numpy as np
import pickle
from typing import List, Dict

class Svm():

    class_code = {'H': '1', 'E': '2', '-': '3'}

    def __init__(self, id_file=False, setype='testset'):
        if id_file:
            with open(id_file) as filein:
                self.dataset = Dict
                self.id_list = filein.read().splitlines()
        else:
            self.dataset = Dict
            self.id_list = []
        
        self.svm_dataset = []
        self.setype = setype

    def encode(self, dataset=False, window=17):
        '''
        Documentation TODO
        '''
        try: dataset != False
        except: 
            print('Method usage: obj.encode(dataset=X)')
            raise SystemExit
        else:
            padding = np.zeros((window//2,20))
            
            for id in self.id_list:

                profile = np.vstack((padding, dataset[id]['profile'], padding))
                dssp = ' '*8 + dataset[id]['dssp'] + ' '*8

                i, j, k = 0, window//2, window
                while k <= len(dssp):
                    row = ''
                    row += self.class_code[dssp[j]]

                    val = 1
                    num_line = i
                    while val <= window * 20:
                        for character in range(20):
                            if profile[num_line][character] == 0.0: 
                                val += 1
                            else:
                                string = ' ' + str(val) + ':' + str(profile[num_line][character])
                                row += string
                                val += 1
                        num_line += 1

                    
                    if self.setype == 'testset' or self.setype == 'blindset':
                        self.svm_dataset.append(row)
                    else:
                        if len(row) > 1:
                            self.svm_dataset.append(row)
                    
                    i, j, k = i+1, j+1, k+1
            
            return(self)

    def fetch_prediction(self, prediction_file):
        try: prediction_file != False
        except: 
            print('Method usage: obj.fetch_prediction(pred_file=X)')
            raise SystemExit
        else:
            with open(prediction_file) as filein:
                prediction = filein.read().splitlines()
                start = 0
                
                for id in self.id_list:
                    lenght_sequence = len(self.dataset[id]['dssp'])
                    end = start + lenght_sequence
                    svm_pred = ''

                    for i in range(start, end):
                        for key,val in self.class_code.items():
                            if val == prediction[i]: svm_pred += key
                    self.dataset[id]['svm_pred'] = svm_pred
                    start = i + 1
        return(self)

    def save(self, path, format) -> object:
        try: path != False; format != False
        except: 
            print('Method usage: obj.save(path=X, format=(dat, pkl))')
            raise SystemExit
        else:
            if format == 'dat':
                with open(path, 'w') as fileout:
                    for i in self.svm_dataset:
                        i += '\n'
                        fileout.write(i)
                    print('%s: done!' %path)

            elif format == 'pkl':
                with open(path, 'wb') as fileout:
                    pickle.dump(self.dataset, fileout)
                    print('%s: done!' %path)

        return(self)

    def load(self, dataset, id_file, encode=False, pkl=False):
        try: dataset != False; id_file != False
        except: 
            print('Method usage: obj.load(path=X, id_file=Y)')
            raise SystemExit
        else:
            with open(id_file) as filein1:
                self.id_list = filein1.read().splitlines()
                if pkl:
                    with open(dataset, 'rb') as filein2:
                        self.dataset = pickle.load(filein2)
                else:
                    self.dataset = dataset

                if encode: self.encode(dataset=self.dataset)
        return(self)

    def fetch_encoded_data(self) -> Dict:
        return(self.svm_dataset)
    
    def fetch_dictionary(self):
        return(self.dataset)