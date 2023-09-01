#takes unselected and selected amino acid sequences, merges into dictionary
#usage: merge.py selected_aminoacids.txt unselected_aminoacids.txt outfile.txt

from collections import defaultdict
from operator import itemgetter
import random
import sys

d1={}
d2={}
d = {}

infile1=sys.argv[1]
with open(infile1) as f_s:
	next(f_s) #this line of code removes the first line of the document, assuming header line
	for line in f_s:
		line=line.strip()
		val_s=line.split('\t')
		seq_s=val_s[0]
		count_s=val_s[1]
		d1[seq_s]=count_s
	#create a dictionary of the sequences and counts in the selected library

infile2=sys.argv[2]		
with open(infile2) as f_un:
	next(f_un) #this line of code removes the first line of the document, assuming header line
	for line in f_un:
		line=line.strip()
		val_un=line.split('\t')
		seq_un=val_un[0]
		count_un=val_un[1]
		d2[seq_un]=count_un
	#create a dictionary of the sequences and counts in the unselected library
		

for key in d2.keys():
	try:
		d[key]=(int(d1[key]),int(d2[key]))
	except KeyError:
		d[key]=(1,int(d2[key]))
#makes new dictionary where the value is a tuple of the selected and unselected counts
#if the sequence is not present in the selected library, sets the count to 1

sort=[]
sort=sorted(d.items(),  key=lambda x:x[1][0], reverse=True)
#sorts by the highest selected count

outfile=sys.argv[3]
with open(outfile, 'w') as f:
	print('Sequence' + '\t' + 'Count' + '\t' + 'Rand',file=f)
with open(outfile, 'a') as f:
	for s in range(len(sort)):
		print(str(sort[s][0]) + "\t" + str(sort[s][1][0]) + "\t" + str(random.randint(110,210)),file=f)
with open(outfile, 'a') as f:
	for s in range(len(sort)):	
		print(str(sort[s][0]) + "\t" + str(sort[s][1][1]) + "\t" + str(random.randint(1,100)),file=f)
#creates an output file with a header line, selected counts with random number between 110 and 210, unselected counts with random number between 1 and 100
