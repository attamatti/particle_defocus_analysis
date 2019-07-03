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

try:
    statedres = float(sys.argv[2])
except:
    sys.exit('USAGE particle_DF_analysis.py <run_data.star> <stated resolution>')

def readfile_get_DF(file):
    '''
    read the run_data.star file with the particles - return list of all defoci
    '''
    defoci = []
    labels,header,data = read_starfile_new(file)
    for i in data:
        DF = (float(i[labels['_rlnDefocusU']])+float(i[labels['_rlnDefocusU']]))/2
        defoci.append(DF)
    return(defoci)

# get yer data
alldefoci = []
for i in sys.argv[1:-1]:
    for j in readfile_get_DF(i):
        alldefoci.append(j)

def calculate_05res(DF):
    '''
    Using y=2.2097x-0.1675 calculated using Henning Stahberg's CTF spreadsheet with 1.0,1.5,2.0,3.0 um DF 
    '''
    return((2.2097*(DF/10000))-0.1675)

pixelsize = 1.065
bincount = {}
bins = [round(pixelsize+(x*pixelsize*0.5),1) for x in range(2,12)]
for i in bins:
    bincount[i] = 0

# calculate 0.05 envelope for all parts check against the bins
srcount = 0 
for i in alldefoci:
    env = calculate_05res(i)
    print (i/10000)
    print(env)
    for j in bins:
        if float(j) > env:
            bincount[j]+=1
    if statedres > env:
        srcount+=1
# plot it
binkeys = list(bincount)
binkeys.sort()
revbinvals = [bincount[x] for x in binkeys]
plt.bar(binkeys,revbinvals,tick_label=binkeys,width=0.2,align='center')
plt.xlabel('Resolution')
plt.ylabel('Effective # particles')
plt.grid(axis='y')
plt.plot(statedres,srcount,marker='o',markersize=10,markeredgecolor='k',c='r')
plt.text(statedres-0.1,srcount,'{0}A\n{1} parts'.format(statedres,srcount))
plt.gca().invert_xaxis()
plt.tight_layout()
plt.savefig('effective_partno.png')
plt.show()