#!/usr/bin/env python

import sys
import matplotlib.pyplot as plt

###---------function: read the star file get the header, labels, and data -------------#######
def read_starfile_new(f):
    inhead = True
    alldata = open(f,'r').readlines()
    labelsdic = {}
    data = []
    header = []
    count = 0
    labcount = 0
    for i in alldata:
        if '_rln' in i:
            labelsdic[i.split()[0]] = labcount
            labcount +=1
        if inhead == True:
            header.append(i.strip("\n"))
            if '_rln' in i and '#' in i and  '_rln' not in alldata[count+1] and '#' not in alldata[count+1]:
                inhead = False
        elif len(i.split())>=1:
            data.append(i.split())
        count +=1
    
    return(labelsdic,header,data)
#---------------------------------------------------------------------------------------------#

def readfile_get_DF(file):
    '''
    read the run_data.star file with the particles - return list of all defoci
    '''
    defoci = []
    labels,header,data = read_starfile_new(file)
    for i in data:
        DF = (float(i[labels['_rlnDefocusU']])+float(i[labels['_rlnDefocusU']]))/20000
        defoci.append(DF)
    return(defoci)
# get yer data
alldefoci = []
for i in sys.argv[1:]:
    for j in readfile_get_DF(i):
        alldefoci.append(j)
bins = [0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0]
n,bins,patches = plt.hist(alldefoci,bins=bins)
plt.xticks(bins)
plt.grid()
plt.savefig('micrograph_DFvals.png')
print(n/max(n))