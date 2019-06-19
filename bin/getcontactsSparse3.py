#!/usr/bin/env python

import os
import sys
from sys import stdout
from Bio.PDB import PDBParser, CaPPBuilder, is_aa
from math import sqrt,pow

# This function takes an input PDB code and returns a structure 
# object.
def pdb_buildstructure(pdbfile):
	pdb_parser = PDBParser(PERMISSIVE=1) 	# The PERMISSIVE instruction allows PDBs presenting errors.
	return pdb_parser.get_structure("name",pdbfile)    # This command gets the structure of the PDB

# This function takes a PDB and returns an object with
# the structural data of chain "chain" of the model "model". 
def pdb_getchain(pdb, chain = "A", model = 0):
	pdb_structure = pdb_buildstructure(pdb)
	return pdb_structure[model][chain]

def norm(atom_1,atom_2):
	return sqrt(pow((atom_1[0]-atom_2[0]),2) + pow((atom_1[1]-atom_2[1]),2) + pow((atom_1[2]-atom_2[2]),2))

# This function outputs a sparse matrix with the residue contacts for
# chain "chain" of the protein "pdb". The variable "dist" is used to
# determine the distance cutoff for residues to be considered in contact.

def pdb_contacts(pdb,chain,dist):	
	i=0
	# Get chain code from 6th letter in pdb name
	pdb_chain = pdb_getchain(pdb,chain)
	ppb=CaPPBuilder()
	# Initialise building of a polypeptide and its sequence
	# If a mutated residue is present in a chain it is classed as a hetatm
	# However, not all hetatms in a chain are part of the sequence. The CaPPBuilder
	# makes sequences by requiring CA-CA distances to be <4.3A. Common hetatms are
	# identified such that an MSE hetatm will be replaced by an M in the sequence
	polypepTot = ppb.build_peptides(pdb_chain, aa_only=False)[0] 
	sequen = polypepTot.get_sequence()

	# Add to the polypeptide
	for polypep_raw in ppb.build_peptides(pdb_chain, aa_only=False)[1:] :
		sequen     += (polypep_raw.get_sequence())
		polypepTot += polypep_raw

	i=0

	# Sometimes the terminal residue in a protein isn't fully resolved
	last_res = polypepTot[-1]
	if last_res.has_id("CA") or last_res.has_id("CB"):
		polypep = polypepTot      # If resolved take whole AA
		file_seq.write(">sequence\n%s\n" %sequen)
#		file_seq.write("%s" %sequen)
	else:
		polypep = polypepTot[:-1] # Otherwise take all but the last AA 	
		file_seq.write(">sequence\n%s\n" %sequen[:-1])		
#		file_seq.write("%s" %sequen[:-1])		

	file_map.write( str(len(polypep)) +"\n" )
#	sys.stderr.write(pdb+'\n')	
	
	for residue1 in polypep:
	# Quite frequently residues do not have resolved CB, in which case use CA
	# If no CA exists, print ERROR. Grep the output if running unsupervised.
			try:	
				if residue1.has_id("CB"): #get_resname() == "GLY":
					c_alpha=residue1["CB"]		
				else:
					c_alpha=residue1["CA"]
			except:
				sys.stdout.write("ERROR")
				raise
			i+=1
			j=0
			for residue2 in polypep:
					try:
						if residue2.has_id("CB"): #get_resname() == "GLY":
							c_alpha2=residue2["CB"]
						else:					
							c_alpha2=residue2["CA"]
					except:
						file_map.write("ERROR")
						raise
					j+=1
					if (norm(c_alpha.get_coord(),c_alpha2.get_coord()) < dist): # 3.5 ):
                                                file_map.write("%d %d\n" %(i-1,j-1))




if __name__ == '__main__':
	argc = len(sys.argv)
	if argc < 2:
                #$QA/bin/getcontactsSparse3.py $TARGET/$DECOY $CHAIN_TARGET 8.0 2>> $TARGET.log

		print "Usage: ./getcontactsSparse3.py <PDBID> Chain 8.0"
	else:
		pdb=sys.argv[1]
		chain=sys.argv[2]
		dist=float(sys.argv[3])
                fasta_extension=".proxy_fasta"
                map_extension=".proxy_map"
                if dist != 8:
                    fasta_extension=fasta_extension+"_"+str(dist)
                    map_extension=map_extension+"_"+str(dist)
                if pdb[-4:] == ".pdb":
                    		file_seq = open(pdb[:-4]+fasta_extension, 'w')
                    		file_map = open(pdb[:-4]+map_extension,'w')
                else:
                                file_seq = open(pdb+fasta_extension, 'w')
                                file_map = open(pdb+map_extension,'w')
                                chain="A"
                pdb_contacts(pdb,chain,dist)

	
