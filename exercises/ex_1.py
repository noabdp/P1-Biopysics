#!/usr/bin/env python
import argparse
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.PDBParser import PDBParser

parser = argparse.ArgumentParser(prog='ex_1', description='Determine the list of pairs of residues whose CA atoms are closer than a given distance')

#add arguments
parser.add_argument('PDB_file' , help='Required PDB file for the program') #PDB file
parser.add_argument('distance',type=float, help='MAX. CA distance in A') #list of distances


# Read command line into args
args = parser.parse_args()


#load structures
pdb_parser =  PDBParser(PERMISSIVE=1)
structure = pdb_parser.get_structure("structure", args.PDB_file) #we load the actual structure from the given pdb

s = []

#get CA's 
for at in structure.get_atoms(): #atoms in our structure that are CA
    if at.id == "CA":
        s.append(at)

#nb search
nbsearch = NeighborSearch(s)

print("CA contacts:")

# Search for contacts
contacts = nbsearch.search_all(args.distance)

# Sort contacts by residue number
contacts.sort(key=lambda pair: (
            pair[0].get_parent().id[1],
            pair[1].get_parent().id[1]
        )
    )

ncontact = 1

for at1, at2 in contacts:

    residue1 = at1.get_parent()
    residue2 = at2.get_parent()

    distance = at1 - at2

    print(
            f"{ncontact:3d}: "
            f"{residue1.get_resname():3s} "
            f"{residue1.id[1]:4d}  --  "
            f"{residue2.get_resname():3s} "
            f"{residue2.id[1]:4d}  "
            f"{distance:6.2f} Å"
        )

    ncontact += 1
