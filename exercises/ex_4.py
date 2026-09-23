#!/usr/bin/env python
import argparse
from Bio.PDB.PDBParser import PDBParser

parser = argparse.ArgumentParser(prog='ex_4', description=' Generate a list of all CA atoms of given residue type with coordinates')

#add arguments
parser.add_argument('PDB_file' , help='Required PDB file for the program') #PDB file
parser.add_argument('res_type',type=str, help='Residue type') 

# Read command line into args
args = parser.parse_args()


#load structures
pdb_parser =  PDBParser(PERMISSIVE=1)
st = pdb_parser.get_structure("structure", args.PDB_file) #we load the actual structure from the given pdb

s = []

for at in st.get_atoms():
    if at.id == "CA":
        parent = at.get_parent()
        if parent.get_resname() == args.res_type.upper(): #get_parent tells us the residue and add upper to have no issues
            s.append(at)


for at in s:
    res = at.get_parent()
    print(f"Residue {res.get_resname()} {res.id[1]} (Chain {res.get_parent().id}): {at.get_coord()}")