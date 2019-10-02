#!/usr/local/bin/python3

from sys import argv 
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


def plot1(dictionary, name, outdir):

    DATAFRAME = pd.DataFrame(frequency(dictionary)).T
    #print(DATAFRAME)
    
    grid_kws = {"height_ratios": (.9, .05), "hspace": .3}
    f, (ax, cbar_ax) = plt.subplots(2, gridspec_kw=grid_kws)
    ax = sns.heatmap(DATAFRAME, ax=ax, cbar_ax=cbar_ax, cbar_kws={"orientation": "horizontal", 'label': 'Residue Frequency'})
    ax.set_title(name)

    plt.savefig(outdir + name +'.png')

def plot2(helix, strand, outdir):

    DF_HELIX = pd.DataFrame(frequency(helix)).T
    DF_STRAND = pd.DataFrame(frequency(strand)).T

    grid_kws = {"height_ratios": (.9, .05), "hspace": .3}
    f, ((ax1, ax2), (cbar_ax1, cbar_ax2)) = plt.subplots(2, 2, gridspec_kw=grid_kws)
    f.suptitle('Residue composition: windows', fontname="DejaVu Sans", fontweight="bold", fontsize=14)

    ax1 = sns.heatmap(DF_HELIX, ax=ax1, cbar_ax=cbar_ax1, cmap="rocket", cbar_kws={"orientation": "horizontal", "label": "Residue frequency"},)
    ax1.set_title("Helix")
    ax2 = sns.heatmap(DF_STRAND, ax=ax2, cbar_ax=cbar_ax2, cmap="rocket", cbar_kws={"orientation": "horizontal", "label": "Residue frequency"},)
    ax2.set_title("Strand")

    plt.savefig(outdir + 'SLW17_Composition.png')

def frequency(dictionary):

    TOT = [0 for i in range(17)]

    for key in dictionary.keys():
        for value in range(17):
            TOT[value] += dictionary[key][value]
    
    for key in dictionary.keys():
        for value in range(17):
            dictionary[key][value] /= TOT[value]
    #print(TOT)
    return(dictionary)


def sliding_window(fasta_file, dssp_file, size):

    RESIDUES = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    STRUCTURE = ['H', 'E']
    HELIX = dict((res,[0 for i in range(17)]) for res in RESIDUES)
    STRAND = dict((res,[0 for i in range(17)]) for res in RESIDUES)
    
    for line in zip(fasta_file, dssp_file):
        if len(line[0]) < 17: continue

        SEQ = line[0].rstrip()
        STR = line[1].rstrip()
        i, j, k = 0, 8, size
        movements = 0

        while k <= len(SEQ):
            if STR[j] not in STRUCTURE: 
                i,j,k = i+1, j+1, k+1
                movements += 1 

            elif STR[j] == 'H':
                for res in range(i, k):
                    RES = SEQ[res]
                    if RES not in RESIDUES: continue
                    HELIX[RES][res - movements] += 1
                i,j,k = i+1, j+1, k+1
                movements += 1

            elif STR[j] == 'E':
                for res in range(i, k):
                    RES = SEQ[res]
                    if RES not in RESIDUES: continue
                    STRAND[RES][res - movements] += 1
                i,j,k = i+1,j+1,k+1
                movements += 1
            
    return(HELIX, STRAND)


if __name__ == '__main__':
    try:
        FASTA_FILE = argv[1] # Results/FASTA_cat.fasta ## FASTA w/o HEADERS FILE!!!
        DSSP_FILE = argv[2] # Results/DSSP_cat.dssp ## DSSP w/o HEADERS FILE!!!
        SIZE = argv[3]
        OUTDIR = argv[4] # Results/
    except:
        print('Program Usage: text.py <FASTA_FILE> <DSSP_FILE> <OUTDIR>')
        raise SystemExit
    else:
        with open(FASTA_FILE) as fasta_file, open(DSSP_FILE) as dssp_file:
            HELIX, STRAND = sliding_window(fasta_file, dssp_file, int(SIZE))
            #df1 = pd.DataFrame(HELIX)
            #print(df1)
        '''
        print()
        frequency(HELIX)
        print()
        frequency(STRAND)
        print()
        '''
        plot2(HELIX, STRAND, OUTDIR)