
# %%
# imports
import MDAnalysis as mda  # for determining whether to construct edges between nodes
import igraph       # for running network analysis on the network at each frame
import numpy as np  # data processing

# %%
# Import data

path = "/home/jnorth/Nextcloud/extracurriculars/publishing/research/hsu/SURE_S2020_fileshare/sure_data/cypA/mdsim_2020_07_10_13_49_07"

u = mda.Universe(path+"md_0_1.gro", path+"md_0_1.xtc")

# convert into a useful format for analysis


# %%
# Write a function that return "1" if two residues are close based on a value "thresh"

def ressAreClose(res1, res2, thresh):
    # if distance between res1 and res2 is less than thresh,
        # return 1
    # else 
    return 0    # return 0

# %%
# Write a function that return "1" if two atoms are close based on a value "thresh"

def atmsAreClose(atom1, atom2, thresh):
    # if distance between atom1 and atom2 is less than thresh,
        # return 1
    # else 
    return 0    # return 0

# %%
# Iterate through a series of frames in an MD trajectory

# Generate a 3D-tensor of dimension (nres x nres x frames) where each element is either 1 or 0 and takes a value based on that timepoint
