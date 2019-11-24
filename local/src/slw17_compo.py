#!/usr/local/bin/python3

import sys, matplotlib.pyplot as plt, matplotlib.ticker as mtick, numpy as np, pandas as pd, seaborn as sns
from Tools import Dataset, Dssp, Fasta, Pssm

residues = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
structures = ['H','E']
pos = ['-8','-7','-6','5','-4','-3','-2','-1','0','1','2','3','4','5','6','7','8']

def slw(dataset, window=17, padding=True, normalize=True):
    tensor = np.zeros((len(structures), window, len(residues)), dtype=np.float64)
    dictionary = {'H': tensor[0], 'E': tensor[1]}
    ss_count = {'H': 0, 'E': 0}
    pad = np.zeros((window//2,20))

    for val in dataset.values():
        fasta = val['fasta']
        dssp = val['dssp']

        try: (len(fasta) + len(dssp)) % 2 == 0
        except:
            print("Different lenght between fasta and dssp: \n%s\n%s" %(fasta, dssp))
            raise SystemExit
        else:
            onehot = np.zeros((len(fasta),len(residues)))
            for res in range(len(fasta)):
                try: residues.index(fasta[res])
                except: continue
                else: onehot[res][residues.index(fasta[res])] = 1.0
           
            if padding: 
                onehot = np.vstack((pad, onehot, pad))
                dssp = ' '*8 + dssp + ' '*8

            i, j, k = 0, window//2, window
            while k <= len(dssp):
                if dssp[j] not in dictionary: 
                    i,j,k = i+1, j+1, k+1
                    continue
                dictionary[dssp[j]] += onehot[i:k]
                ss_count[dssp[j]] += 1
                i,j,k = i+1, j+1, k+1
    
    if normalize:
        normalizer = 0
        for i in range(len(structures)):
            normalizer += ss_count[structures[i]]
        
        for i in range(len(structures)):
            dictionary[structures[i]] /= normalizer
            ss_count[structures[i]] /= normalizer
    
    return(dictionary)


def slw_heatmap(dictionary, outdir=False, RGB='Greys_r'):
    color = RGB
    df_helix = pd.DataFrame(dictionary['H']).T
    df_strand = pd.DataFrame(dictionary['E']).T

    df_helix.columns = pos
    df_helix.index = residues
    df_strand.columns = pos
    df_strand.index = residues

    grid_kws = {"height_ratios": (.9, .05), "hspace": .3}
    f, ((ax1, ax2), (cbar_ax1, cbar_ax2)) = plt.subplots(2, 2, gridspec_kw=grid_kws)
    f.suptitle('Residue composition: windows17', fontname="Times New Roman", fontweight="bold", fontsize=9)

    ax1 = sns.heatmap(df_helix, ax=ax1, cbar_ax=cbar_ax1, cmap=color, cbar_kws={"orientation": "horizontal", "label": "Residue frequency"},)
    ax1.set_title("Helix", fontname="Times New Roman", fontweight="bold", fontsize=9)
    ax2 = sns.heatmap(df_strand, ax=ax2, cbar_ax=cbar_ax2, cmap=color, cbar_kws={"orientation": "horizontal", "label": "Residue frequency"},)
    ax2.set_title("Strand", fontname="Times New Roman", fontweight="bold", fontsize=9)

    if outdir:
        plt.savefig(outdir + 'slw17_compo.png', box_inches='tight')
    else:
        plt.show()


if __name__ == '__main__':
    try:
        setype = sys.argv[1]
    except:
        print('Program Usage: res_compo <setype (trainingset, blindset)>')
        raise SystemExit
    else:
        wd = 'dataset/' + setype + '/'
        if setype == 'trainingset':
            data_id = wd + 'ts2.id'
        elif setype == 'blindset':
            data_id = wd + 'bs.id'

        # Build the dataset from scratch
        fasta = Fasta(data_id, setype=setype).parse(fastatype='singleline').fetch_dict()
        dssp = Dssp(data_id, setype=setype, raw_file=False).parse().fetch_dict()
        dataset = Dataset(data_id, setype=setype).build(fasta=fasta, dssp=dssp).fetch_dict()

        dictionary = slw(dataset)
        slw_heatmap(dictionary)