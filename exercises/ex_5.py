#!/usr/bin/env python
import argparse
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.PDBParser import PDBParser

parser = argparse.ArgumentParser(prog='ex_5', description="Generate a list of backbone connectivity (peptide bonds between C and N atoms)")

#add arguments
parser.add_argument('PDB_file' , help='Required PDB file for the program') #PDB file
parser.add_argument("cutoff", type=float,default =2.5, help='Cut-off distance for peptide bonds in Å') 


# Read command line into args
args = parser.parse_args()


#load structures
pdb_parser =  PDBParser(PERMISSIVE=1)
structure = pdb_parser.get_structure("structure", args.PDB_file) #we load the actual structure from the given pdb

s = []
ordinary = {"C","N"}


for at in structure.get_atoms():
    if at.id in ordinary:
        s.append(at)

#nb search
nbsearch = NeighborSearch(s)
# Search for contacts
contacts = nbsearch.search_all(args.cutoff)

p_bonds = []

for at1, at2 in contacts:
    res1 = at1.get_parent() #get residues
    res2 = at2.get_parent() #get other residue

    # ignore same res
    if res1 == res2:
        continue

    if {at1.id, at2.id} == ordinary:
        distance = at1 - at2 
        p_bonds.append((at1, at2, distance))

print(f"Backbone Peptide Bonds: ")
ncontact = 1

for at1, at2, dist in p_bonds:
    res1 = at1.get_parent()
    res2 = at2.get_parent()

    print(f"{ncontact:}: "
        f"{res1.get_resname():3s} ({at1.id}) -- "
        f"{res2.get_resname():3s} ({at2.id})"
        f"{dist:6.2f} Å")
    ncontact += 1