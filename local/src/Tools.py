#!/usr/local/bin/python3
import numpy as np, pickle, sys
from typing import List, Dict
from sys import argv

class Dataset():
    '''
    The class Dataset is intended to be updated as soon as some new formatted files 
    will be introduced in my bioinformatics pipelines. For the time being, sequence profiles
    coming from psiblast runs and dssp outputs are parsed both starting from raw and pre-pruned 
    files, then datasets are built by stacking information of single objects (by the use of other classes)
    in order to deal with the whole DATABASE when training models in a Machine Learning framemwork.
    '''
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

    def parse(self, normalize=False):
        '''
        Documentation TODO
        '''
        onehot = []
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
                            if res not in self.residues: continue
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
                                np.save('./dataset/' + self.setype + '/psiblast/bin/' + id, self.dataset[id])
                        else:
                            self.dataset[id] = array
                            np.save('./dataset/' + self.setype + '/psiblast/bin/' + id, self.dataset[id])

            else:
                self.dataset[id] = np.load(path_to_profile + 'bin/' + id + '.npy', allow_pickle=True)
                if normalize:
                    self.dataset[id] /= 100
                if np.sum(self.dataset[id]) == 0:
                    onehot.append(id)
                    with open('./dataset/' + self.setype  + '/fasta/' + id + '.fasta') as fasta_file:
                                sequence = fasta_file.read().splitlines()[1]
                                for res in range(len(sequence)):
                                    try: self.residues.index(sequence[res])
                                    except: continue
                                    else: self.dataset[id][res][self.residues.index(sequence[res])] = 1.0
        if len(onehot) > 0:
            print('These profiles have been one-hot encoded!\n', onehot)

        return(self)

    def save(self):
        for key in self.dataset.keys():
            np.save('./dataset/' + self.setype + '/psiblast/bin/' + key, self.dataset[key])

    def fetch_dict(self) -> Dict:
        return(self.dataset)


class Dssp():

    def __init__(self, id_file, setype, raw_file=False, pdb=False):
        with open(id_file) as filein:
            id_list = filein.read().splitlines()
        
        self.residues = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
        self.raw_file = raw_file
        self.setype = setype
        self.pdb = pdb
        self.dataset = dict((id, False) for id in id_list)

    def parse(self):
        try:
            self.raw_file == False
        except:
            print('Method usage: obj.parse(chain=X)')
            raise SystemExit
        else:
            dictionary = {  'H': 'H', 'G': 'H', 'I': 'H', 
                            'E': 'E', 'B': 'E', 
                            'C': '-', 'S': '-', 'T': '-', ' ': '-'}
            unknown = []
            for key in self.dataset.keys():
                if self.pdb: 
                    id = key.split('_')[0]
                    ch = key.split('_')[1]
                if self.raw_file: 
                    try: 
                        dssp_file = open('dataset/' +  self.setype + '/dssp/raw/' + id + '.dssp')
                    except: 
                        unknown.append(key)
                        continue
                    else:
                        ss = ''
                        flag = 0
                        for row in dssp_file:
                            if row.find('  #  RESIDUE') == 0:
                                flag = 1
                                continue
                            if flag == 1:
                                if row[13] == '!': continue
                                if row[11] == ch:
                                    ss += dictionary[row[16]]

                        self.dataset[key] = ss
                else:
                    try: dssp_file = open('dataset/' + self.setype + '/dssp/ss/' + key + '.dssp')
                    except: 
                        unknown.append(key)
                        continue
                    else: self.dataset[key] = dssp_file.read().splitlines()[1]
            
            if len(unknown) > 0:
                print('These DSSP have not been found!\n', unknown)

        return(self)
    
    def save(self):
        for key in self.dataset.keys():
            with open('dataset/' + self.setype + '/dssp/ss/' + key + '.dssp','w') as fileout:
                fileout.write('>' + key + '\n' + self.dataset[key] + '\n') 

    def fetch_dict(self) -> Dict:
        return(self.dataset)

class Fasta():
    def __init__(self, id_file, setype):
        with open(id_file) as filein:
            id_list = filein.read().splitlines()
        
        self.residues = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
        self.setype = setype
        self.dataset = dict((id, False) for id in id_list)

    def parse(self, fastatype=False):
        try:
            fastatype == False
        except:
            print('Method usage: obj.parse(fastatype=(singleline, multiline))')
            raise SystemExit
        else:
            for id in self.dataset:
                path_to_fasta = './dataset/' + self.setype + '/fasta/'
                if fastatype == 'singleline':
                    with open(path_to_fasta + id + '.fasta') as fasta_file:
                        self.dataset[id] = fasta_file.read().splitlines()[1]
        return(self)

    def fetch_dict(self) -> Dict:
        return(self.dataset)

if __name__ == '__main__':
    id_list = sys.argv[1]
    setype = sys.argv[2]
    prof = Pssm(id_list, setype=setype, raw_file=True).parse().save()
    # dict_prof = prof.fetch_dict()

    #dssp = Dssp(id_list, setype=setype, raw_file=True, pdb=True).parse().save()
    #ls dict_dssp = dssp.fetch_dict()
    # data = Dataset(id_list, setype='testset').build(profile=dict_prof, dssp=dict_dssp)
    # dataset = data.fetch_dict()
    # for key in dataset:
    #     print(key, '\n', dataset[key]['profile'], '\n')#, dataset[key]['dssp'])

    # data.save(path='./dataset/cv/fold5/', name='cv4.pkl')
