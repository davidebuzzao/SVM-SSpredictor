#!/usr/local/bin/python3
import sys, matplotlib.pyplot as plt, numpy as np, pandas as pd
from matplotlib import cm

def taxa_piechart(dictionary, outdir=False, other=False, RGB='Greys_r'):
    # Pie chart, where the slices will be ordered and plotted counter-clockwise:
        labels = []
        sizes = []
        color = RGB


        if other == True:
            for key,val in dictionary.items():
                if val[1] < 5:
                    labels.append(key)
                    sizes.append(val[1])
        else:
            other = 0
            for key,val in dictionary.items():
                if val[1] >= 5:
                    labels.append(key)
                    sizes.append(val[1])
                else:
                    other += val[1]
            labels.append('Other')
            sizes.append(other)

        fig1, ax1 = plt.subplots()
        theme = plt.get_cmap(color)
        ax1.set_prop_cycle("color", [theme(1. * i / len(sizes))
                                    for i in range(len(sizes))])
        if other == True:
            ax1.pie(sizes,
                    labeldistance=1.2,
                    startangle=70,
                    wedgeprops={"edgecolor":"k", 'linewidth': 1, 'linestyle': 'solid', 'antialiased': True})
        else:
            explosion = np.zeros(len(labels))
            explosion[labels.index('Other')] = 0.05
            ax1.pie(sizes,
                    explode=explosion,
                    labeldistance=1.2,
                    startangle=70,
                    wedgeprops={"edgecolor":"k", 'linewidth': 1, 'linestyle': 'solid', 'antialiased': True})
        
        plt.legend( loc='upper left',
                    labels=['%s, %1.1f%%' % (
                            l, s) for l, s in zip(labels, sizes)],
                    prop={'size': 9},
                    bbox_to_anchor=(0.0, 1),
                    bbox_transform=fig1.transFigure)

        #draw circle
        centre_circle = plt.Circle((0,0),0.50,fc='white', edgecolor='black', linewidth=1)
        fig = plt.gcf()
        fig.gca().add_artist(centre_circle)

        # Equal aspect ratio ensures that pie is drawn as a circle
        ax1.axis('equal')
        plt.title('Taxa Composition', fontname="Times New Roman", fontweight="bold")
        if outdir:
            plt.savefig(outdir + 'taxa_compo.png', box_inches='tight')
        else:
            plt.show()


if __name__ == '__main__':
    try:
        table = sys.argv[1]
    except:
        print('Program Usage: taxa_compo <domains_table>')
        raise SystemExit
    else:
        dictionary = {}
        with open(table) as filein:
            for line in filein:
                line = line.rstrip().split()
                dictionary[line[0]] = dictionary.get(line[0],[int(line[1]),float(line[2])])

    df = pd.DataFrame(dictionary).T
    df.columns = ['Frequency','Percentage']

    taxa_piechart(dictionary)
    taxa_piechart(dictionary)