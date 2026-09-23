#!/usr/bin/env python
import argparse
from Bio.PDB.PDBParser import PDBParser

parser = argparse.ArgumentParser(prog='ex_7', description="Distances between all atom pairs of two given residues")

#add arguments
parser.add_argument('PDB_file' , help='Required PDB file for the program') #PDB file
parser.add_argument("res_1", type=int) 
parser.add_argument("res_2", type=int) 


# Read command line into args
args = parser.parse_args()


#load structures
pdb_parser =  PDBParser(PERMISSIVE=1)
st = pdb_parser.get_structure("structure", args.PDB_file) #we load the actual structure from the given pdb

s = []
g_res = {args.res_1,args.res_2}

for at in st.get_atoms():
    res = at.get_parent()
    if res.id[1] in g_res:
        s.append(at)

#separate them
res1_atoms = []
res2_atoms = []

for at in s:
    res = at.get_parent()

    if res.id[1] == args.res_1:
        res1_atoms.append(at)
    elif res.id[1] == args.res_2:
        res2_atoms.append(at)

#distances

for at1 in res1_atoms:
    for at2 in res2_atoms:
        distance = at1 - at2

        print(f"{at1.id} -- {at2.id} "
            f"{distance:6.2f} Å"
        )