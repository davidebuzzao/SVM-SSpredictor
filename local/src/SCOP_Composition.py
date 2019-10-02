#!/usr/local/bin/python3

import sys
import matplotlib.pyplot as plt

def pretty_dictionary(D):
    print()
    for key,val in D.items():
        print(key, " = ", val)
    print()

def unknownDOM(unknown, OUT_DIR):
    with open(OUT_DIR, 'w') as OUTFILE:
        for element in unknown:
            OUTFILE.write(element+'\n')

## Build the legend SCOP-classes dictionary 
def SCOP_dict(scop_des_file):
    LEGEND = {}
    with open(scop_des_file) as FILE_IN:
        for line in FILE_IN:
            LINE = line.split()
            ID = LINE[0]
            CLASS = ' '.join(LINE[1:])
            if len(ID) == 1:
                LEGEND[ID] = LEGEND.get(ID, [CLASS,0])
    return(LEGEND)

## Computation of SCOP-class composition 
def SCOP_compo(LEGEND, dom_file, scop_class_file):
    with open(dom_file) as DOM_FILE, open(scop_class_file) as SCOP_FILE:
        DICT_ID = dict([(domain.rstrip(),'') for domain in DOM_FILE])
        
        for line in SCOP_FILE:
            ID = str(line.split()[0].rstrip())
            if ID in DICT_ID:
                DICT_ID[ID] = line.split()[1].split('.')[0]
        
        UNKNOWN = []
        TOT = 0
        for key,val in DICT_ID.items():
            if val in LEGEND:
                LEGEND[val][1] += 1
                TOT += 1
            else:
                UNKNOWN.append(key)

    return(DICT_ID, LEGEND, TOT, UNKNOWN)
        

if __name__ == '__main__':
    try:
        DOM_FILE = sys.argv[1] ## jpred4.list.txt
        LEGEND_FILE = sys.argv[2] ## legend.txt (retrieved by cleaning dir.des.scope.2.06-stable.txt)
        SCOP_CLASS_FILE = sys.argv[3] ## dir.cla.scope.2.06-stable.txt
        #OUT_DIR = sys.argv[4] ## ./Results/
    except:
        print('Program Usage: script.py <DOM_FILE> <CLEAN_SCOP_DOM_FILE> <SCOP_CLASS_FILE> <OUT_DIR>')
        raise SystemExit
    else:
        LEGEND1 = SCOP_dict(LEGEND_FILE)
        ID, LEGEND2, TOT, UNKNOWN = SCOP_compo(LEGEND1, DOM_FILE, SCOP_CLASS_FILE)
       
        pretty_dictionary(LEGEND2)
        #unknownDOM(UNKNOWN, OUT_DIR+'unknown_domains.txt')

        # Pie chart, where the slices will be ordered and plotted counter-clockwise:
        labels = []
        sizes = []
        for key,val in LEGEND2.items():
            if val[1] != 0:
                labels.append(val[0])
                sizes.append(round((val[1]/TOT)*100,2))

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
        plt.title('SCOP_composition', fontname="Times New Roman", fontweight="bold")
        #plt.show()
        #plt.savefig(OUT_DIR + 'SCOP_composition', box_inches='tight')