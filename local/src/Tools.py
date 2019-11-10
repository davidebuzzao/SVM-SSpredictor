#!/usr/local/bin/python3
from sys import argv
import numpy as np
import pickle
from typing import List, Dict

class Dataset():

    info_stored: List
    dataset: Dict
    
    def __init__(self, id_file=False, setype=False):
        self.setype = setype
        if id_file:
            with open(id_file) as filein:
                id_list = filein.read().splitlines()
        
            info_stored = [ 'profile',
                            'dssp',
                            'fasta',
                            'gor_pred',
                            'svm_pred',
                            'ss_compo',
                            'r_compo']

            self.dataset = dict((id, dict((info, False) for info in info_stored)) for id in id_list)
        else:
            self.dataset = None

    def build(self, profile=False, dssp=False, fasta=False):
        if profile:
            for key,val in profile.items():
                self.dataset[key]['profile'] = val
        if dssp:
            for key,val in dssp.items():
                self.dataset[key]['dssp'] = val
        if fasta:
            for key,val in fasta.items():
                self.dataset[key]['fasta'] = val

        return(self)

    def save(self, path=False, name=False) -> object:
        try: path != False; name != False
        except: 
            print('Method usage: obj.save(path=X, name=Y)')
            raise SystemExit
        else:
            with open(path + name, 'wb') as fileout:
                pickle.dump(self.dataset, fileout)
        return(self)

    def load(self, path=False):
        try: path != False
        except: 
            print('Method usage: obj.load(path=X)')
            raise SystemExit
        else:
            with open(path, 'rb') as filein:
                self.dataset = pickle.load(filein)
        return(self)

    def fetch_dict(self):
        return(self.dataset)
    
    def __iter__(self):
        return(self)

    def __len__(self):
        return len(self.dataset)


class Pssm():

    def __init__(self, id_file, setype, raw_file=False):
        with open(id_file) as filein:
            id_list = filein.read().splitlines()

        self.residues = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
        self.raw_file = raw_file
        self.setype = setype
        self.dataset = dict((id, False) for id in id_list)

    def parse(self):
        '''
        Documentation TODO
        '''
        for id in self.dataset:
            init_profile = []
            path_to_profile = './dataset/' +  self.setype + '/psiblast/'
            if self.raw_file: 
                try:
                    pssm_file = open(path_to_profile + 'pssm/' + id + '.pssm')
                except:
                    with open('./dataset/' + self.setype + '/fasta/'  + id + '.fasta') as fasta_file:
                        sequence = fasta_file.read().splitlines()[1]
                        for res in sequence:
                            prof = np.zeros(20)
                            try: self.residues.index(sequence[res])
                            except: continue
                            else: 
                                prof[self.residues.index(res)] = 1.0
                                init_profile.append(prof)
                    self.dataset[id] = np.array(init_profile, dtype=np.float64)
                    np.save('./dataset/' + self.setype + '/psiblast/bin/' + id, self.dataset[id])
                else:
                    for row in pssm_file:
                        row = row.rstrip().split()
                        try: row[0]
                        except: continue
                        else:
                            if row[0].isdigit():
                                init_profile.append(row[22:-2])
                    array = np.array(init_profile, dtype=np.float64)
                    
                    if self.setype == 'testset' or self.setype == 'blindset':
                        if np.sum(array) == 0:
                            with open('./dataset/' + self.setype  + '/fasta/' + id + '.fasta') as fasta_file:
                                sequence = fasta_file.read().splitlines()[1]
                                for res in range(len(sequence)):
                                    try: self.residues.index(sequence[res])
                                    except: continue
                                    else: array[res][self.residues.index(sequence[res])] = 1.0
                                self.dataset[id] = array
                                np.save('./dataset/' + self.setype + '/psiblast/bin/' + id, array)
                        else:
                            self.dataset[id] = array/100
                            np.save('./dataset/' + self.setype + '/psiblast/bin/' + id, array)

            else:   
                self.dataset[id] = np.load(path_to_profile + 'bin/' + id + '.npy', allow_pickle=True)
                if np.sum(self.dataset[id]) == 0:
                    with open('./dataset/' + self.setype  + '/fasta/' + id + '.fasta') as fasta_file:
                                sequence = fasta_file.read().splitlines()[1]
                                for res in range(len(sequence)):
                                    try: self.residues.index(sequence[res])
                                    except: continue
                                    else: self.dataset[id][res][self.residues.index(sequence[res])] = 1.0
        return(self)

    def fetch_dict(self) -> Dict:
        return(self.dataset)


class Dssp():

    def __init__(self, id_file, setype, raw_file=False):
        with open(id_file) as filein:
            id_list = filein.read().splitlines()
        
        self.residues = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
        self.raw_file = raw_file
        self.setype = setype
        self.dataset = dict((id, False) for id in id_list)

    def parse(self, chain=False):
        try:
            self.raw_file == False
        except:
            print('Method usage: obj.parse(chain=X)')
            raise SystemExit
        else:
            dictionary = {  'H': 'H', 'G': 'H', 'I': 'H', 
                            'E': 'E', 'B': 'E', 
                            'C': '-', 'S': '-', 'T': '-', ' ': '-'}
            for id in self.dataset:
                path_to_dssp = './dataset/' +  self.setype + '/dssp'
                if self.raw_file:
                    with open(path_to_dssp + '/raw/' + id + '.dssp') as dssp_file:
                        flag = 0
                        for row in dssp_file:
                            if row.find('  #  RESIDUE') == 0:
                                flag = 1
                                continue
                            if flag == 1:
                                if row[13] == '!': continue
                                if row[11] == self.chain:
                                    self.secondary_structure += dictionary[row[16]]
                else:
                    with open(path_to_dssp + '/ss/' + id + '.dssp') as dssp_file:
                        self.dataset[id] = dssp_file.read().splitlines()[1]

        return(self)

    def fetch_dict(self) -> Dict:
        return(self.dataset)

if __name__ == '__main__':
    id_list = argv[1]
    prof = Pssm(id_list, setype='testset', raw_file=True).parse()
    dict_prof = prof.fetch_dict()

    dssp = Dssp(id_list, setype='testset', raw_file=False).parse()
    dict_dssp = dssp.fetch_dict()

    data = Dataset(id_list, setype='testset').build(profile=dict_prof, dssp=dict_dssp)
    dataset = data.fetch_dict()
    # for key in dataset:
    #     print(key, '\n', dataset[key]['profile'], '\n')#, dataset[key]['dssp'])

    data.save(path='./dataset/cv/fold5/', name='cv4.pkl')
