#!/usr/local/bin/python3
import sys
import numpy as np

if __name__ == '__main__':
    filein = sys.argv[1]
    list_id = ''
    
    with open(filein) as f:
        for id in f:
            id = id.rstrip()
            
            path = ''
            try: prof = np.load('./psiblast/bin/' + id + '.npy')
            except: pass
            else:
                if np.sum(prof) == 0:
                    list_id += id + '\n'
    print(list_id)