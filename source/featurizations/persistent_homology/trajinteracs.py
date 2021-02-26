
'''
Determines how interactions between proteins and ligands evolve over time during the course of a trajectory. May use getcontacts! (great package for measuring interactions)
'''


'''
IMPORTS
'''

import Bio


'''
FUNCTIONS
'''

# distance cutoff
def distCutoffPlot(atom1: Bio.PDB.Atom.Atom, atom2:Bio.PDB.Atom.Atom, thresh: float) -> bool:
    '''
    If the distance is less than the threshold, return 0. If greater, return 1
    '''

    # euclidean norm

    return




# plot distances between atoms of interest


# compute the contact matrix for a given interaction at a timestep



