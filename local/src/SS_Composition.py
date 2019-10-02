#!/usr/local/bin/python3

import sys 
import matplotlib.pyplot as plt
import argparse

def ss_compo(dssp_file):
    ss = ['H', 'E', '-']
    dictionary={'H':0, 'E':0, '-':0}
    for line in dssp_file:
        line = str(line.rstrip())
        for ch in line:
            dictionary[ch] += 1

    return (dictionary['H'], dictionary['E'], dictionary['-'])

if __name__ == '__main__':
    # try:
    #     dssp = sys.argv[1]
    #     #outdir = sys.argv[2]
    # except:
    #     print('Program Usage: text.py <DSSP_FILE> <OUT_DIR>')
    #     raise SystemExit
    # else:
    parser = argparse.ArgumentParser(description='Compute Confusion Matrix performance and FOM')
    parser.add_argument('dssp_file', help = 'This is the dssp file with concatenating BlindSet secondary structure')
    parser.add_argument('--dest', help = 'This is the relative/absolute path to your outuput folder')
    args = parser.parse_args()
    outdir = args.dest 
    dssp = args.dssp_file
    with open(dssp) as filein:
        H, E, C = ss_compo(filein)
        TOT = H + E + C
        print(TOT)

    # Pie chart, where the slices will be ordered and plotted counter-clockwise:
    labels = 'Helix', 'Strand', 'Coil'
    sizes = [H/TOT, E/TOT, C/TOT]

    #colors
    colors = ['#ff9999','#66b3ff','#99ff99']
    
    fig1, ax1 = plt.subplots()
    ax1.pie(sizes, colors = colors, labels=labels, autopct='%1.1f%%', startangle=90)
    
    #draw circle
    centre_circle = plt.Circle((0,0),0.70,fc='white')
    fig = plt.gcf()
    fig.gca().add_artist(centre_circle)

    # Equal aspect ratio ensures that pie is drawn as a circle
    ax1.axis('equal')
    plt.title('SS composition', fontname="Times New Roman", fontweight="bold")
    #plt.show()
    plt.savefig(outdir + '/SS_composition.png', box_inches='tight')
