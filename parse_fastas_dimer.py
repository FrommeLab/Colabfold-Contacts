"""
Takes an input fasta file and writes out input files for local colabfold
The bait sequence should be named >bait. Use this version if the bait is a dimer

"""

import argparse
from Bio import SeqIO
from Bio import Seq 


parser = argparse.ArgumentParser()
parser.add_argument('--i', type=str, required=True, help='input file')               
args = parser.parse_args()

input_file = open(args.i)

def load_dict():
	my_dict = SeqIO.to_dict(SeqIO.parse(input_file, "fasta"))
	return my_dict




f = load_dict()
x = f.keys()
bait = f['bait']
#inserted extra bait copy
for i in x:
	t = f[i]
	combined = "id,sequence\n"+t.id.replace('.','_')+','+bait.seq+':'+bait.seq+':'+t.seq
	file = open(t.id.replace('.','_') + ".csv", "w")
	file.write(str(combined))
	file.close()
	
