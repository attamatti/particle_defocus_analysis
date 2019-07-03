#!/usr/bin/env python

# new version with nice graphs!

import sys
import matplotlib.pyplot as plt 
import numpy as np
###---------function: read the star file get the header, labels, and data -------------#######
def read_starfile(f):
    alldata = open(f,'r').readlines()
    labelsdic = {}
    data = []
    header = []
    for i in alldata:
        if '#' in i:
            labelsdic[i.split('#')[0]] = int(i.split('#')[1])-1
        if len(i.split()) > 3:
            data.append(i.split())
        if len(i.split()) < 3:
            header.append(i.strip("\n"))
    return(labelsdic,header,data)
#---------------------------------------------------------------------------------------------#

if len(sys.argv) != 2:
    sys.exit("USAGE: starfile_get_defocus_range.py <star file>")

(labels,header,data) = read_starfile(sys.argv[1])
defoci = []
for i in data:
    defocus = (float(i[labels['_rlnDefocusU ']]) + float(i[labels['_rlnDefocusV ']]))/20000
    defoci.append(defocus)

bins = [0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0]
n,bins,patches = plt.hist(defoci,bins=bins)
plt.xticks(bins)
plt.grid()   
plt.savefig('particle_DFrange.png')
print '''
-------
Results
-------
{0} items

Min defocus:  {1}
Max defocus:  {2}
mean defocus: {3}
std           {4}
'''.format(len(defoci),min(defoci),max(defoci),np.mean(defoci),np.std(defoci))
print(n/max(n))
