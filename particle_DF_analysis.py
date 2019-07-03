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
    cutoff = sys.argv[3]
    pixelsize = 1.065
except:
    sys.exit('USAGE particle_DF_analysis.py <run_data.star> <stated resolution> <cutoff 0.05 or 0.01>')
if cutoff not in ['0.05','0.01']:
    sys.exit('USAGE particle_DF_analysis.py <run_data.star> <stated resolution> <cutoff 0.05 or 0.01>\nERROR: invalid cutoff value - use 0.01 or 0.05')

#-----------------------------------------------------------------------#
def calculate_0501res(DF,val):
    '''
    Using y = 2.2097x-0.1675 for 0.5
    and
    y = 1.3885x - 0.0173 for 0.01
    calculated using Henning Stahberg's CTF spreadsheet with 0.5,1.5,2.0,3.0 um DF 
    '''
    if val == '0.05':
        return((2.2097*(DF/10000))-0.1675)
    if val == '0.01':
        return((1.3885*(DF/10000))-00.0173)
#-----------------------------------------------------------------------#

# make the dict of bins
bincount = {}
bins = [round(pixelsize+(x*pixelsize*0.5),1) for x in range(2,12)]
for i in bins:
    bincount[i] = 0
goodparts = []

# calculate 0.05 envelope for all parts check against the bins
labels,header,data = read_starfile_new(sys.argv[1])
for i in data:
    defocus = (float(i[labels['_rlnDefocusU']])+float(i[labels['_rlnDefocusU']]))/2
    env = calculate_0501res(defocus,cutoff)
    for j in bins:
        if float(j) > env:
            bincount[j]+=1
    if statedres > env:
        goodparts.append(i)
        
# plot it
binkeys = list(bincount)
binkeys.sort()
revbinvals = [bincount[x] for x in binkeys]
plt.bar(binkeys,revbinvals,tick_label=binkeys,width=0.2,align='center')
plt.xlabel('Resolution')
plt.ylabel('Effective # particles')
plt.grid(axis='y')
plt.plot(statedres,len(goodparts),marker='o',markersize=10,markeredgecolor='k',c='r')
plt.text(statedres-0.1,len(goodparts),'{0}A\n{1} parts'.format(statedres,len(goodparts)))
plt.gca().invert_xaxis()
plt.title('cutoff = {0}'.format(cutoff))
plt.tight_layout()
plt.savefig('effective_partno.png')
plt.show()

# write a starfile with selected particles
starout = open('CTFe_filt.star','w')
for i in header:
    starout.write('{0}\n'.format(i))
for i in goodparts:
    starout.write('{0}\n'.format('    '.join(i)))