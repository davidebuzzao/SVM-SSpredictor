#!/usr/local/bin/python3
from sys import argv
import numpy as np

class Database:
    '''
    The class Database is intended to be updated as soon as some new formatted files 
    will be introduced in my bioinformatics pipelines. For the time being, sequence profiles
    coming from psiblast runs and dssp outputs are parsed both starting from raw and pre-pruned 
    files, then datasets are built by stacking information of single objects (by the use of other classes)
    in order to deal with the whole DATABASE when training models in a Machine Learning framemwork.
    '''
    dataset: list
    datatype: str
    raw_file: bool
    window: int
    one_hot: bool
    
    def __init__(self, datatype, setype=False, psiblast_raw_file=False, dssp_raw_file=False, window=17, one_hot=False):
        try:
            datatype != False
        except:
            print('Method Usage: obj.build_dataset(id_list, datatype = (psiblast, dssp, svm)')
            raise SystemExit
        else:
            self.residues = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
            self.dataset = []
            self.datatype = datatype
            self.psiblast_raw_file = psiblast_raw_file
            self.dssp_raw_file = dssp_raw_file
            self.window = window
            self.one_hot = one_hot
            self.setype = setype

    def __len__(self):
        return len(self.dataset) 
    
    def __append__(self, item):
        return self.dataset.append(item)

    def build_dataset(self, id_list):

        if self.datatype == 'psiblast':
            for id in id_list:
                if self.psiblast_raw_file:
                    if self.one_hot:
                        id = Pssm('./psiblast/pssm/' + id + '.pssm', './fasta/' + id + '.fasta', id=id, raw_file=self.psiblast_raw_file, setype=self.setype)
                    else:
                        id = Pssm('./psiblast/pssm/' + id + '.pssm', id=id, raw_file=self.psiblast_raw_file)
                else:
                    id = Pssm('./psiblast/bin/' + id + '.npy')

                self.__append__(id.parser())
            return self.dataset

        elif self.datatype == 'dssp':
            for id in id_list:

                if self.dssp_raw_file:
                    ch = id.rstrip().split('_')[1]
                    id = id.rstrip().split('_')[0]
                    id = Dssp(path_dssp='./dssp/raw/' + id + '.dssp', chain=ch, dssp_raw_file=True)
                else:
                    id = Dssp(path_dssp='./dssp/ss/' + id + '.dssp', dssp_raw_file=False)
                self.__append__(id.parser())
            return self.dataset

        elif self.datatype == 'svm':
            svm = Svm(id_list, psiblast_raw_file=self.psiblast_raw_file, dssp_raw_file=self.dssp_raw_file, window=self.window, setype=self.setype)
            svm_dataset = svm.encode()  
            return(svm_dataset)

    def __iter__(self):
        return(self)

############################################################
############################################################

class Pssm(Database):

    path_profile: str

    def __init__(self, path_profile, path_fasta=False, id=False, raw_file=False, setype=False):
        super().__init__(self)
        
        self.id = id
        self.path_profile = path_profile
        self.path_fasta = path_fasta
        self.profile = np.empty(0)
        self.raw_file = raw_file
        self.setype = setype

    def parser(self):
        init_profile = []
        if self.raw_file:
            try: 
                pssm_file = open(self.path_profile)
            except:
                # print(self.id)
                with open(self.path_fasta) as fasta_file:
                    sequence = fasta_file.read().splitlines()[1]
                    for res in sequence:
                        prof1 = np.zeros(20)
                        prof1[self.residues.index(res)] = 1.0
                        init_profile.append(prof1)
                self.profile = np.array(init_profile, dtype=np.float64)
                np.save('./psiblast/bin/' + self.id + '.npy', self.profile)
            else:
                for row in pssm_file:
                    row = row.rstrip().split()
                    try: row[0]
                    except: pass
                    else:
                        if row[0].isdigit():
                            init_profile.append(row[22:-2])
                self.profile = np.array(init_profile, dtype=np.float64)
                self.profile/100
                np.save('./psiblast/bin/' + self.id + '.npy', self.profile)
                
                if self.setype == 'test':
                    
                    with open(self.path_fasta) as fasta_file:
                        sequence = fasta_file.read().splitlines()[1]

                        for res in range(len(sequence)):
                            if np.sum(self.profile[res]) == 0:
                                try: self.residues.index(sequence[res])
                                except: pass
                                else:
                                    # print(self.id)
                                    self.profile[res][self.residues.index(sequence[res])] = 1.0
        else:   
            self.profile = np.load(self.path_profile, allow_pickle=True)
            if np.sum(self.profile) == 0:
                with open(self.path_fasta) as fasta_file:
                        sequence = fasta_file.read().splitlines()[1]

                        for res in range(len(sequence)):
                            if np.sum(self.profile[res]) == 0:
                                try: self.residues.index(sequence[res])
                                except: pass
                                else:
                                    # print(self.id)
                                    self.profile[res][self.residues.index(sequence[res])] = 1.0
        
        return(self.profile)
    
    def __len__(self):
        return len(self.dataset)

############################################################
############################################################

class Dssp(Database):
    
    path_dssp: str
    chain: str
    secondary_structure: str

    def __init__(self, path_dssp, chain=False, dssp_raw_file=False):
        super().__init__(self)
        self.path_dssp = path_dssp
        self.chain = chain
        self.secondary_structure = ''
        self.dssp_raw_file = dssp_raw_file

    def parser(self):
        dictionary = {
            'H': 'H', 'G': 'H', 'I': 'H', 
            'E': 'E', 'B': 'E', 
            'C': '-', 'S': '-', 'T': '-', ' ': '-'
            }
        
        if self.dssp_raw_file:
            with open(self.path_dssp) as dssp_file:
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
            with open(self.path_dssp) as dssp_file:
                self.secondary_structure = dssp_file.read().splitlines()[1]

        return(self.secondary_structure)

    def ss_composition(self):
        ss = ['H', 'E', '-']
        frequencies = {'H':0, 'E':0, '-':0}

        for sublist in self.dataset:
            for ch in sublist:
                frequencies[ch] += 1

        return frequencies
    
    def __len__(self):
        return len(self.dataset)

############################################################
############################################################

class Svm(Database):

    def __init__(self, id_list, psiblast_raw_file, dssp_raw_file, window, setype):
        super().__init__(self)
        
        self.id_list = id_list
        self.psiblast_raw_file = psiblast_raw_file
        self.dssp_raw_file = dssp_raw_file
        self.window = window
        self.setype = setype

        self.svm_data = []
        self.padding = np.zeros((self.window//2,20))
        self.dictionary = {'H': '1', 'E': '2', '-': '3'}

    
    def encode(self):
        '''
        The main idea:
        - build profiles and dssps dataset
        1. iterate for every couple in the datasets
            - padding operation both for prof and dssp is performed
            2. iterate a sliding window N*20 'til the lenght of the padded dssp
                3. iterate a val in a range of N*20
                    4. iterate 20 times (every column of a profile's row) for every line in the window[i:k]
                        - if the column is zero skip
                        - else save the column in the specified format (val:column)
        '''
        profiles = Database(datatype='psiblast', setype='test', psiblast_raw_file=self.psiblast_raw_file, window=17, one_hot=True)
        profiles.build_dataset(self.id_list)
        dataset_profile = profiles.dataset

        dssps = Database(datatype='dssp', dssp_raw_file=False)
        dssps.build_dataset(self.id_list)
        dataset_dssp = dssps.dataset

        for item in range(len(dataset_dssp)):
            profile = np.vstack((self.padding, dataset_profile[item], self.padding))
            dssp = ' '*8 + dataset_dssp[item] + ' '*8

            i, j, k = 0, self.window//2, self.window
            while k <= len(dssp):
                row = ''
                row += self.dictionary[dssp[j]]

                val = 1
                num_line = i
                while val <= self.window * 20:
                    for character in range(20):
                        if profile[num_line][character] == 0.0: 
                            val += 1
                        else:
                            string = ' ' + str(val) + ':' + str(profile[num_line][character])
                            row += string
                            val += 1
                    num_line += 1
                if len(row) > 1:
                    self.svm_data.append(row)
                i, j, k = i+1, j+1, k+1
        
        return self.svm_data

############################################################
############################################################

if __name__ == '__main__':
    id_input = argv[1]

    id_list = []
    with open(id_input) as filein:
        for id in filein:
            id_list.append(id.rstrip())
    
    prof = Database(datatype='svm', setype='test', psiblast_raw_file=True, dssp_raw_file=False, window=17, one_hot=True)
    prof_dataset = prof.build_dataset(id_list)

    for i in prof_dataset:
        print(i)
        # for j in i:
            # print(j)

    # svm = Database(datatype='svm', raw_file=False, window=17, one_hot=True)
    # svm_dataset = svm.build_dataset(id_input)
    # for i in svm_dataset:
    #     print(i)