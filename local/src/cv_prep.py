#!/usr/local/bin/python3

'''
The program takes in input a randomnly shuffled Y ids list (via sort -R program), 
the path of a folder and the number X of folds to perform Cross-Validation analysis.
The program outputs X id lists in the folder specified and every file has Y/X ids. 
'''

import sys

if __name__ == '__main__':
    try:
        input_id = sys.argv[1]
        directory = sys.argv[2]
        folds = int(sys.argv[3])
    except:
        print('Program Usage: python3 cv_prep.py <Input random id list> <Output Folder> <Number of FOLDS>')
        raise SystemExit
    else:
        with open(input_id) as filein:
            name = filein.readlines()
            block = int(len(name)/folds)
            output_id = [directory + 'cv' + str(i) + '.id' for i in range(1,folds+1)]

            i = 0
            for j in range(folds):
                with open(output_id[j], 'w') as fileout:
                    for i in range(i,i+block):
                        fileout.write(name[i])