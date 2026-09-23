#!/usr/bin/env python
import argparse
from Bio.PDB.PDBParser import PDBParser

parser = argparse.ArgumentParser(prog='ex_2', description=' Generate a list of all atoms for a given residue number')

#add arguments
parser.add_argument('PDB_file' , help='Required PDB file for the program') #PDB file
parser.add_argument('residue',type=int, help='Residue number (Including Chain if applicable)') 
parser.add_argument("--chain", default="A", help="Chain identifier")

# Read command line into args
args = parser.parse_args()


#load structures
pdb_parser =  PDBParser(PERMISSIVE=1)
st = pdb_parser.get_structure("structure", args.PDB_file) #we load the actual structure from the given pdb

s = []

for at in st.get_atoms():
    res = at.get_parent()          # parent residue object
    chain = res.get_parent()       # parent chain object
    
    # Filter
    if res.id[1] == args.residue and chain.id == args.chain:
        s.append(at)

if s:
    res_name = s[0].get_parent().get_resname()
    print(f"Residue {res_name} {args.residue} (Chain {args.chain}):")
    for atom in s:
        print(f"  Atom: {atom.get_name():<4} Coordinates: {atom.get_coord()}")
else:
    print(f"No atoms found for Residue {args.residue} on Chain {args.chain}.")