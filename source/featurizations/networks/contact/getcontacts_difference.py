
# %%
# imports


import igraph           # igraph to calculate network parameters
import numpy as np      # numpy to manipulate matrices
import Bio.PDB.PDBParser as PDBParser
from Bio.PDB.Polypeptide import PPBuilder

# %%


# get a maximal RIN for the protein
# res x res
# each res is of the format RES_###

p = PDBParser()
structure = p.get_structure("1a2k", "data/JLN_1_21_1C/batch-download-structures-1608267263488/5q0i.pdb")[0]["A"]
ppb=PPBuilder()

seq = str()

for pp in ppb.build_peptides(structure):
    #print(pp.get_sequence())
    seq += pp.get_sequence()

#seq = pp.get_sequence().__str__()

# %%

# get unbound structure residue interaction network adjacency matrix


# fill in a copy of the maximal RIN with obtained residue interactions
len(seq)

# get bound structure RIN adj matrix


# fill in a copy of the maximal RIN with obtained residue interactions


# Subtract the two matrices


# delete the rows & columns that have only zeroes


# visualize the non-zero differences

