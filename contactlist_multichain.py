#!/usr/bin/env python
"""
List contacts <12 angstroms between two chains
Author: Ryan Feathers jrf296
Date: 11/22/2022
Chains can be specified with --id1 and --id2 
"""
import argparse
import os.path
import pandas as pd
import csv
import Bio.PDB
import numpy
from Bio.PDB import PDBParser, PDBIO
def calc_residue_dist(residue_one, residue_two) :
    """Returns the C-alpha distance between two residues"""
    diff_vector  = residue_one["CA"].coord - residue_two["CA"].coord
    return numpy.sqrt(numpy.sum(diff_vector * diff_vector))

def calc_dist_matrix(chain_one, chain_two) :
    """Returns a matrix of C-alpha distances between two chains"""
    answer = numpy.zeros((len(chain_one), len(chain_two)), float)
    for row, residue_one in enumerate(chain_one) :
        for col, residue_two in enumerate(chain_two) :
            answer[row, col] = calc_residue_dist(residue_one, residue_two)
    return answer

parser = argparse.ArgumentParser()
parser.add_argument('--o', type=str, required=False,  nargs='?', const=1, 
	default='scores.csv' , help='output file name')                 
parser.add_argument('--id1', type=int, required=False,  nargs='?', const=1, 
	default=0 , help='chain 1 index')
parser.add_argument('--id2', type=int, required=False,  nargs='?', const=1, 
	default=1 , help='chain 2 index')	 
args = parser.parse_args()

#makes a list of files for iterative use
path ='./'
files = os.listdir(path)
pdb_l = []
scores = []
dm = []

for i in files:
	if 'pdb' in str(i):
		pdb_l.append(i) 
for i in pdb_l:
	pdb_code = i
	pdb_filename = i


	chain_ids = []
	pdb = PDBParser().get_structure(pdb_code, pdb_filename)

	for chain in pdb.get_chains():
		chain_ids.append(chain.get_id())
	
	if len(chain_ids) > 2:
		print("There are more than 2 chains in the model use index args")	

	structure = Bio.PDB.PDBParser().get_structure(pdb_code, pdb_filename)
	model = structure[0]

	dist_matrix = calc_dist_matrix(model[str(chain_ids[args.id1])], model[str(chain_ids[args.id2])])
	contact_map = dist_matrix < 12

	df = pd.DataFrame(dist_matrix)

	contacts = []
	cutoff = 12

	for i in range(df.shape[0]): #iterate over rows
	    for j in range(df.shape[1]): #iterate over columns
	        value = df.at[i, j] #get cell value
	        if value <= cutoff:
	        	contacts.append([i+1,j+1,value])

	#this is not perfect, reversing drop dups from B to A gives different results
	contacts_matrix = pd.DataFrame(contacts)
	if contacts_matrix.empty:
		scores.append([str(pdb_filename),' '+ 'contacts = '+'0'])
		print(str(pdb_filename),' '+ 'contacts = '+'0')
	else:
		contacts_matrix.columns = ['ChainA','ChainB','Distance']
		contacts_matrix = contacts_matrix.sort_values('Distance',ignore_index=True)
		contacts_matrix = contacts_matrix.drop_duplicates('ChainA', keep='first', ignore_index=True)
		contacts_matrix = contacts_matrix.drop_duplicates('ChainB', keep='first', ignore_index=True)
		scores.append([str(pdb_filename),' '+ 'contacts = '+str(len(contacts_matrix.index))])
		print(str(pdb_filename),' '+ 'contacts = '+str(len(contacts_matrix.index)))
		dm.append([str(pdb_filename), pd.DataFrame(dist_matrix)])
#still need to scores.sort based on len(scores[1])

 
	# opening the csv file in 'w+' mode
file = open(args.o, 'w+', newline ='')
 
	# writing the data into the file
with file:   
	write = csv.writer(file)
	write.writerows(scores)
