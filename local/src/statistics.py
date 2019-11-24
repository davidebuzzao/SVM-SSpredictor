#!/usr/local/bin/python3

from Tools import Dataset, Dssp, Fasta, Pssm
import sys, matplotlib.pyplot as plt, numpy as np, pandas as pd, seaborn as sns
from res_compo import rSS_compo, rSS_stat, rSS_barplot, ss_piechart
from slw17_compo import slw, slw_heatmap
from scop_compo import scop_legend, scop_compo, scop_piechart
from taxa_compo import taxa_piechart

if __name__ == '__main__':
    try:
        setype = sys.argv[1]
        color = sys.argv[2] 
    except:
        print('Program Usage: res_compo <setype (trainingset, blindset), color(default "Greys_r")>')
        raise SystemExit
    else:
        wd = 'dataset/' + setype + '/'
        db = 'dataset/db/'
        if setype == 'trainingset':
            data_id = wd + 'ts.id'
            table = wd + 'ts.taxa'
        elif setype == 'blindset':
            data_id = wd + 'bs_mary.id'
            table = wd + 'bs_mary.taxa'
        elif setype == 'onehot':
            data_id = wd + 'onehot.id'
            table = wd + 'onehot.taxa'
        
        # Build the dataset from scratch
        fasta = Fasta(data_id, setype=setype).parse(fastatype='singleline').fetch_dict()
        dssp = Dssp(data_id, setype=setype, raw_file=False).parse().fetch_dict()
        dataset = Dataset(data_id, setype=setype).build(fasta=fasta, dssp=dssp).fetch_dict()

        ## rSS statistics
        res_dict, str_dict = rSS_compo(dataset)
        dataframe1 = rSS_stat(res_dict,str_dict)
        rSS_barplot(dataframe1, setype, RGB=color)
        ss_piechart(str_dict, setype, RGB=color)

        ## slw17 statistics 
        dictionary1 = slw(dataset)
        slw_heatmap(dictionary1, RGB=color)

        ## SCOP statistics
        scop_legend_file = db + 'dir.des.scope.2.06-stable.txt' ## legend.txt (retrieved by cleaning dir.des.scope.2.06-stable.txt)
        scop_class_file = db + 'dir.cla.scope.2.06-stable.txt' ## dir.cla.scope.2.06-stable.txt
        legend = scop_legend(scop_legend_file)
        id, dictionary2, tot, unknown = scop_compo(data_id, legend, scop_class_file)
        scop_piechart(dictionary2, tot, RGB=color)

        ## taxa statistics 
        dictionary3 = {}
        with open(table) as filein:
            for line in filein:
                line = line.rstrip().split()
                dictionary3[line[0]] = dictionary3.get(line[0],[int(line[1]),float(line[2])])
        df = pd.DataFrame(dictionary3).T
        df.columns = ['Frequency','Percentage']

        taxa_piechart(dictionary3, RGB=color)
        taxa_piechart(dictionary3, other=True, RGB=color)