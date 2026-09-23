#!/usr/bin/env python
import argparse
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.PDBParser import PDBParser

parser = argparse.ArgumentParser(prog='ex_3', description="Determine all possible hydrogen bonds")

#add arguments
parser.add_argument('PDB_file' , help='Required PDB file for the program') #PDB file
parser.add_argument('distance',type=float, default=3.5, help='Max distance in Å, default: 3.5 Å') #list of distances


# Read command line into args
args = parser.parse_args()


#load structures
pdb_parser =  PDBParser(PERMISSIVE=1)
structure = pdb_parser.get_structure("structure", args.PDB_file) #we load the actual structure from the given pdb

s = []
polar = {"0", "N", "S"}

#get all polar atoms
for at in structure.get_atoms(): 
    if at.element in polar:
        s.append(at)

#nb search
nbsearch = NeighborSearch(s)

print(f"Possible Hydrogen Bonds:")

# Search for contacts
contacts = nbsearch.search_all(args.distance)

# sort contacts by residue number
contacts.sort(
    key=lambda pair: (
        pair[0].get_parent().get_parent().id,
        pair[0].get_parent().id[1],
        pair[1].get_parent().get_parent().id,
        pair[1].get_parent().id[1],
    )
)

ncontact = 1

for at1, at2 in contacts:
    res1 = at1.get_parent()
    res2 = at2.get_parent()

    # skip same residue contacts 
    if res1 == res2:
        continue

    chain1 = res1.get_parent().id
    chain2 = res2.get_parent().id
    distance = at1 - at2

    print(f"{ncontact:}: "
        f"({at1.get_name().strip()}) -- "
        f"({at2.get_name().strip()})  "
        f"{distance:6.2f} Å")

    ncontact += 1
