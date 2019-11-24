#!/usr/local/bin/python3
import sys, matplotlib.pyplot as plt, numpy as np
from matplotlib import cm

def pretty_dictionary(D):
    print()
    for key,val in D.items():
        print(key, " = ", val)
    print()

def unknownDOM(unknown, OUT_DIR):
    with open(OUT_DIR, 'w') as outfile:
        for element in unknown:
            outfile.write(element+'\n')

## Build the legend SCOP-classes dictionary 
def scop_legend(scop_des_file):
    legend = {}
    with open(scop_des_file) as filein:
        for line in filein:
            line = line.rstrip().split()
            if line[0] == '#': continue
            id = line[2]
            if len(id) == 1:
                cla = ' '.join(line[4:])
                legend[id] = legend.get(id, [cla,0])
    return(legend)

## Computation of SCOP-class composition 
def scop_compo(dom_file, legend, scop_class_file):
    with open(dom_file) as filein1, open(scop_class_file) as filein2:
        dict_id = dict([(domain.rstrip(),'') for domain in filein1])

        for line in filein2:
            line = line.rstrip().split()
            id = line[0]
            if id in dict_id:
                dict_id[id] = line[3].split('.')[0]

        unknown = []
        tot = 0
        for key,val in dict_id.items():
            if val in legend:
                legend[val][1] += 1
                tot += 1
            else:
                unknown.append(key)
    return(dict_id, legend, tot, unknown)

def scop_piechart(dictionary, tot, outdir=False, RGB='Greys_r'):
    # Pie chart, where the slices will be ordered and plotted counter-clockwise:
        color = RGB
        labels = []
        sizes = []
        for key,val in dictionary.items():
            if val[1] != 0:
                labels.append(val[0])
                sizes.append(round((val[1]/tot)*100,2))

        fig1, ax1 = plt.subplots()
        theme = plt.get_cmap(color)
        ax1.set_prop_cycle("color", [theme(1. * i / len(sizes))
                                    for i in range(len(sizes))])
        
        ax1.pie(sizes,
                #explode=(0.05,0.05,0.05,0.05,0.05,0.05,0.05),
                #labels=labels,
                labeldistance=1.2,
                #autopct='%1.1f%%',
                startangle=70,
                wedgeprops={"edgecolor":"k", 'linewidth': 1, 'linestyle': 'solid', 'antialiased': True})
        
        plt.legend( loc='upper left',
                    labels=['%s, %1.1f%%' % (
                            l, s) for l, s in zip(labels, sizes)],
                    prop={'size': 9},
                    bbox_to_anchor=(0.0, 1),
                    bbox_transform=fig1.transFigure)

        #draw circle
        centre_circle = plt.Circle((0,0),0.50,fc='white', edgecolor='black')
        fig = plt.gcf()
        fig.gca().add_artist(centre_circle)

        # Equal aspect ratio ensures that pie is drawn as a circle
        ax1.axis('equal')
        plt.title('SCOP Composition', fontname="Times New Roman", fontweight="bold")
        if outdir:
            plt.savefig(outdir + 'scop_compo.png', box_inches='tight')
        else:
            plt.show()


if __name__ == '__main__':
    try:
        dom_file = sys.argv[1] ## jpred4.list.txt
        legend_file = sys.argv[2] ## legend.txt (retrieved by cleaning dir.des.scope.2.06-stable.txt)
        scop_class_file = sys.argv[3] ## dir.cla.scope.2.06-stable.txt
    except:
        print('Program Usage: script.py <DOM_FILE> <CLEAN_SCOP_DOM_FILE> <SCOP_CLASS_FILE> <OUT_DIR>')
        raise SystemExit
    else:
        legend = scop_legend(legend_file)
        id, dictionary, tot, unknown = scop_compo(dom_file, legend, scop_class_file)
       
        pretty_dictionary(dictionary)
        scop_piechart(dictionary, tot)