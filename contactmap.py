#!/usr/bin/env python
"""
Generate a CSV of the number of contacts <12A between 2 chains
Author: Ryan Feathers jrf296
Date: 08/25/2022
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
args = parser.parse_args()

#makes a list of files for iterative use
path ='./'
files = os.listdir(path)
pdb_l = []
scores = []

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
		print("There are more than 2 chains in the model")	



	structure = Bio.PDB.PDBParser().get_structure(pdb_code, pdb_filename)
	model = structure[0]

	dist_matrix = calc_dist_matrix(model[str(chain_ids[0])], model[str(chain_ids[1])])
	#contact_map is not actually used
  contact_map = dist_matrix < 12

	df = pd.DataFrame(dist_matrix)

	contacts = []
	cutoff = 12

	for i in range(df.shape[0]): #iterate over rows
	    for j in range(df.shape[1]): #iterate over columns
	        value = df.at[i, j] #get cell value
	        if value <= cutoff:
	        	contacts.append([i+1,j+1,value])

	#this is not perfect, reversing drop dups from B to A gives different results?
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


 
	# opening the csv file in 'w+' mode
file = open(args.o, 'w+', newline ='')
 
	# writing the data into the file
with file:   
	write = csv.writer(file)
	write.writerows(scores)




