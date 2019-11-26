#!/usr/local/bin/python3
import numpy as np, sys, pandas as pd, matplotlib.pyplot as plt, seaborn as sns, matplotlib.ticker as mtick
from Tools import Dataset, Dssp

def fragments(seq, ss):
    fragments = []
    val = 0
    while val < len(seq):
        tmp_list = []
        while seq[val] == ss and val < len(seq) - 1:
            tmp_list.append(val)
            val += 1
        if len(tmp_list) > 0:
            fragments.append(set(tmp_list))
        val += 1

    return(fragments)

if __name__ == '__main__':
    try:
        data_id = sys.argv[1] 
        setype = sys.argv[2]
    except:
        print('Program Usage: gor_training <cv_train.txt>')
        raise SystemExit
    else:
        secondary_structure = ['-', 'H', 'E']
        dssp = Dssp(data_id, setype=setype, raw_file=False).parse()
        dict_dssp = dssp.fetch_dict()

        dict_lenght = {'H':{}, 'E':{}, '-':{}}
        ss_frag_freq = {'H':0, 'E':0, '-':0}
        for key in dict_dssp:
            for ss in secondary_structure:
                if ss in dict_dssp[key]:
                    frag_list = fragments(dict_dssp[key], ss)
                    for frag in frag_list:
                        dict_lenght[ss][len(frag)] = dict_lenght[ss].get(len(frag), 0)
                        dict_lenght[ss][len(frag)] += 1
                        ss_frag_freq[ss] += 1
                else: continue
        print(ss_frag_freq)
        for ch in secondary_structure:
            for key,val in dict_lenght.items():
                for seg in val:
                    val[seg] /= ss_frag_freq[ch]
                    
        df = pd.DataFrame(dict_lenght)
        for ch in secondary_structure:
            df[ch] /= np.sum(df[ch]) * 1/100
            df[ch] = round(df[ch],2)
        df = df.plot(kind='bar', rot=0, colormap='plasma', edgecolor='black')

        df.set(xlabel="Segment Lenght", ylabel="Segment Lenghts Frequency", title='SS Segment Lenghts Composition')
        df.yaxis.set_major_formatter(mtick.PercentFormatter())
        plt.show()