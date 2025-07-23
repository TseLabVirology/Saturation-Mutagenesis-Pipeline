#takes unselected and selected amino acid sequences, calculates percent identity to wt sequence, merges into dictionary
#usage: merge.py selected_aminoacids.txt unselected_aminoacids.txt outfile.txt

from collections import defaultdict
from operator import itemgetter
import random
import sys
from thefuzz import fuzz

infile1=sys.argv[1]
infile2=sys.argv[2]
outfile=sys.argv[3]


wt="DVEGCH"

d0={}
dU={}
dU0 = {}


with open(infile1) as f_s0:
    next(f_s0) #this line of code removes the first line of the document, assuming header line
    for line in f_s0:
        line=line.strip()
        val_s0=line.split('\t')
        seq_s0=val_s0[0]
        m_s0 = fuzz.ratio(wt,seq_s0)
        count_s0=val_s0[1]
        d0[seq_s0]=(count_s0,m_s0)
    #create a dictionary of the sequences and counts in the selected library



with open(infile2) as f_un:
    next(f_un) #this line of code removes the first line of the document, assuming header line
    for line in f_un:
        line=line.strip()
        val_un=line.split('\t')
        seq_un=val_un[0]
        m_un = fuzz.ratio(wt,seq_un)
        count_un=val_un[1]
        dU[seq_un]=(count_un,m_un)


for key in dU.keys():
    try:
        dU0[key]=(int(d0[key][0]),int(dU[key][0]),dU[key][1])
    except KeyError:
        dU0[key]=(0.1,int(dU[key][0]),dU[key][1])

for key in d0.keys():
    try:
        dU0[key]=(int(d0[key][0]),int(dU[key][0]),dU[key][1])
    except KeyError:
        dU0[key]=(int(d0[key][0]),0.1,d0[key][1])
#makes new dictionary where the value is a tuple of the selected and unselected counts
#if the sequence is not present in the selected library, sets the count to 1


sort0=[]
sort0=sorted(dU0.items(),  key=lambda x:x[1][0], reverse=True)
#sorts by the highest selected count

with open(outfile, 'w') as f:
    print('Sequence' + '\t' + 'Count' + '\t' + 'Fuzzy' + '\t' + 'Rand',file=f)
with open(outfile, 'a') as f:
    for s in range(len(sort0)):    
        print(str(sort0[s][0]) + "\t" + str(sort0[s][1][1]) + "\t" + str(sort0[s][1][2]) + "\t" + str(random.randint(5,95)),file=f)
        
with open(outfile, 'a') as f:
    for s in range(len(sort0)):
        print(str(sort0[s][0]) + "\t" + str(sort0[s][1][0]) + "\t" + str(sort0[s][1][2]) + "\t" + str(random.randint(105,195)),file=f)

        
#creates an output file with a header line, selected counts with random number between 110 and 210, unselected counts with random number between 1 and 100

