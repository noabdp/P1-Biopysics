#!/usr/bin/env python
import argparse
import sys
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.PDBParser import PDBParser

parser = argparse.ArgumentParser(prog="ex_6",
    description="Find S-S contacts (disulfide bonds) between CYS residues")

#add arguments
parser.add_argument('PDB_file' , help='Required PDB file for the program') #PDB file
parser.add_argument("cutoff", type=float,  nargs="?", default = 2.5, help='Cut-off distance for disulfide bonds in Å') #more than 1.9 to allow variability


# Read command line into args
args = parser.parse_args()


#load structures
pdb_parser =  PDBParser(PERMISSIVE=1)
structure = pdb_parser.get_structure("structure", args.PDB_file) #we load the actual structure from the given pdb

s = []

#get sulfur atoms from CYS 
for at in structure.get_atoms(): 
    res = at.get_parent()
    if res.get_resname() == "CYS" and at.id == "SG" :
        s.append(at)

if not s:
    print(f"No CYS residues with S atoms found in {args.PDB_file}.")
    sys.exit()

#nb search
nbsearch = NeighborSearch(s)
# Search for contacts
contacts = nbsearch.search_all(args.cutoff)

print(f"Disulfide bonds:")

ncontact = 1
for at1, at2 in contacts:
    res1 = at1.get_parent()
    res2 = at2.get_parent()

    # ignore same res
    if res1 == res2:
        continue

    distance = at1 - at2

    print(
        f"{ncontact:}: "
        f"{res1.get_resname():3s} ({at1.id}) -- "
        f"{res2.get_resname():3s} ({at2.id}) "
        f"{distance:6.2f} Å"
    )

    ncontact += 1

if ncontact == 1:
    print("No disulfide bonds found within the given cutoff distance.")