#!/usr/local/bin/python3

import sys 
import matplotlib.pyplot as plt
import numpy as np

def pretty_dictionary(D):
    for key,val in D.items():
        print(key, "\t=\t", val)

def pretty_matrix(M):
    for sublist in M:
        print(str(sublist[0]) + ':' + str(round(sublist[2], 2)))

def parsingFASTA_UNIPROT(file_in):
    DICT_OS = {}
    for line in file_in:
        if line[0] != '>': continue
        OS = ' '.join(line.split('OS=')[1].split(' ')[0:2])
        
        DICT_OS[OS] = DICT_OS.get(OS, 0)
        DICT_OS[OS] += 1
    
    return DICT_OS

if __name__ == '__main__':
    try:
        FASTA_FILE = sys.argv[1]
        OUT_DIR = sys.argv[2] ## ./Results/
    except:
        print('Program Usage: text.py <FASTA_FILE> <OUT_DIR>')
        raise SystemExit
    else:
        with open(FASTA_FILE) as FILE_IN:
            ORGANISMS = parsingFASTA_UNIPROT(FILE_IN)

        TOT = 0
        for key,val in ORGANISMS.items():
            TOT += val
            '''
            OS_LIST = []
            for key,val in ORGANISMS.items():
                OS_LIST.append([key, val, (val/TOT)*100])
            '''
            #pretty_dictionary(ORGANISMS)
            
        # Pie chart, where the slices will be ordered and plotted counter-clockwise:
        labels = []
        sizes = []
        OTHERS = 0

        for key,val in ORGANISMS.items():
            PERCENTAGE = (val/TOT)*100
            if PERCENTAGE > 2.0:
                labels.append(key)
                sizes.append(PERCENTAGE)
            else:
                OTHERS += PERCENTAGE
        
        labels.append('Others')
        sizes.append(OTHERS)

        ## colors
        colors = ['#ff9999','#66b3ff','#99ff99']
        
        fig1, ax1 = plt.subplots()
        ax1.pie(sizes, colors = colors, labels=labels, autopct='%1.1f%%', startangle=90)
        
        #draw circle
        centre_circle = plt.Circle((0,0),0.70,fc='white')
        fig = plt.gcf()
        fig.gca().add_artist(centre_circle)

        # Equal aspect ratio ensures that pie is drawn as a circle
        ax1.axis('equal')
        plt.title('TAXA_composition', fontname="Times New Roman", fontweight="bold")
        #plt.show()
        plt.savefig(OUT_DIR + 'TAXA_composition', box_inches='tight')