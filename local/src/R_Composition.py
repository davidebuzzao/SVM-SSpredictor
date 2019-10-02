#!/usr/local/bin/python3

import sys 
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def pretty_dictionary(D):
    for key,val in D.items():
        print(key, " = ", val)
    print()

def pretty_matrix(M):
    for sublist in M:
        print(sublist)

### Total composition 1st function
def TOTr_compo(fasta_file):
    residues = {}
    for line in fasta_file:
        line = str(line.rstrip())
        for ch in line:
            residues[ch] = residues.get(ch, 0) ## when trying to use get option and there's no recognised key, 
            residues[ch] += 1                  ## a new key is created with a 0 value initialized
    return residues


########################################
### Total and Fractional composition ###
########################################

## Fastest way to build the dictionary 

def rSS_DICT(fasta_file, dssp_file): ## files in input are given in sorted order

    DOMAINS = {}
    ID = ''

    for line in zip(fasta_file, dssp_file): ## each line becomes a tuple

        if line[0].startswith('>'):
            ID = line[0].rstrip()                  ## Being the files in input sorted in a lexicographycal order, 
            DOMAINS[ID] = DOMAINS.get(ID, ['','']) ## the tuple with domains contain 2 equal values: one of the two will become ID
        else:
            DOMAINS[ID] = [line[0].rstrip(),line[1].rstrip()] ## For every ID a list value with sequence and structure is saved
           
    return DOMAINS

## Very long way to build the dictionary
def rSS_DICT2(fasta_file, dssp_file):

    DOMAINS = {}
    header = ''
    sequence = ''
    structure = ''

    for line in fasta_file:
        line = line.rstrip()
        if header and sequence:
            if line[0] == '>':
                
                DOMAINS[header][0] = sequence
                header = ''
                sequence = ''
            
        if line[0] == '>':
            ID = line
            DOMAINS[ID] = DOMAINS.get(ID, ['',''])
            header = line

        else:
            if header:
                sequence = line
        
    if header and sequence:
        DOMAINS[header][0] = sequence
        header = ''
    
    for line in dssp_file:
        line = line.rstrip()
        if header and structure:
            if line[0] == '>':
                DOMAINS[header][1] = structure
                header = ''
                structure = ''
            
        if line[0] == '>':
            ID = line
            DOMAINS[ID] = DOMAINS.get(ID, ['',''])
            header = line

        else:
            if header:
                structure = line
        
    if header and structure:
        DOMAINS[header][1] = structure

    return DOMAINS
    
## Computation of Residue composition 
def rSS_compo(domains):

    RESIDUES = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    #STRUCTURES = ['B', 'C', 'E', 'G', 'I', 'H', 'S', 'T' ]
    STRUCTURES = ['C', 'E', 'H']
   
    RES_DICT = dict([(res, [0.0, 0.0, 0.0, 0.0]) for res in RESIDUES]) ## Dictionary format = {res: [num_res, H, E, C]}
    STRUCT_DICT = dict([(struct, 0.0) for struct in STRUCTURES]) ## Dictionary format = {SS: num_SS}

    for key,val in DOMAINS.items():
        sequence = val[0]
        structure = val[1]

        for character in range(len(sequence)):
            
            if sequence[character] not in RES_DICT: continue ## just for working with the traditional 20 aa
            RES_DICT[sequence[character]][0] += 1

            if structure[character] == 'H':
                RES_DICT[sequence[character]][1] += 1
                STRUCT_DICT['H'] += 1
            
            elif structure[character] == 'E':
                RES_DICT[sequence[character]][2] += 1
                STRUCT_DICT['E'] += 1
            
            elif structure[character] == '-':
                RES_DICT[sequence[character]][3] += 1
                STRUCT_DICT['C'] += 1

    return(RES_DICT, STRUCT_DICT)
        
if __name__ == '__main__':
    try:
        FASTA_FILE = sys.argv[1]
        DSSP_FILE = sys.argv[2]
        OUT_DIR = sys.argv[3]
    except:
        print('Program Usage: text.py <FASTA_FILE> <DSSP_FILE> <OUT_DIR>')
        raise SystemExit
    else:
        with open(FASTA_FILE) as FASTA_IN, open(DSSP_FILE) as DSSP_IN:
            #RESIDUES = TOTr_compo(FASTA_IN)
            DOMAINS = rSS_DICT(FASTA_IN, DSSP_IN)
            #pretty_dictionary(DOMAINS)
            RES_DICT, STRUCT_DICT = rSS_compo(DOMAINS)
            
            print('\nTOT Residue composition, H, E, C:')
            pretty_dictionary(RES_DICT)

            print('TOT SS composition (H, E, C):')
            pretty_dictionary(STRUCT_DICT)

        ### Overall and SS residue statistics
        TOT = 0
        for key,val in STRUCT_DICT.items():
            TOT += val
        
        PERCENTAGES = {} 
        for key,val in RES_DICT.items():
            PERCENTAGES[key] = PERCENTAGES.get(key, [round((val[0]/TOT)*100,2), round((val[1]/STRUCT_DICT['H'])*100,2), round((val[2]/STRUCT_DICT['E'])*100,2), round((val[3]/STRUCT_DICT['C'])*100,2)])
        DATAFRAME = pd.DataFrame(PERCENTAGES).T
        DATAFRAME.rename(index = {0:'Overall', 1:'Helix', 2:'Strand', 3:'Coil'}, inplace = True) # rename rows indexes

        DATAFRAME.plot.bar(rot = 0) # use a method for dataframe
        plt.title('Residue_Composition', fontname="Times New Roman", fontweight="bold")
        plt.savefig(OUT_DIR + 'Residue_Composition', box_inches='tight')
        #print(DATAFRAME)
        #print('Percentages (TOT, H, E, C):')
        #pretty_dictionary(PERCENTAGES)

        ###########################################################
        ###             PROPENSITY definition                   ###
        ### P(residue, sec_struct) / P(residue) * P(sec_struct) ###
        ###########################################################
        ## Here not computing probabilities
        PROPENSITY = {} 
        for key,val in RES_DICT.items():
            FR_TOT = val[1] + val[2] + val[3]
            PROPENSITY[key] = PROPENSITY.get(key, [round((val[1]/FR_TOT)*100,2), round((val[2]/FR_TOT)*100,2), round((val[3]/FR_TOT)*100,2)])
        
        #print('Propensity Scale Percentage (H, E, C):')
        #pretty_dictionary(PROPENSITY)