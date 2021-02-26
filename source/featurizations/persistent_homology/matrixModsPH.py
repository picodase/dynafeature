

'''
This package includes modifications to matrices to be used in persistent homology computations. For example, custom distance matrices will be implemented to enable selection of specific interaction types.

'''

'''
IMPORTS
'''
import numpy as np
import Bio.PDB


'''
FUNCTIONS
'''

def isBonded(atom1: Bio.PDB.Atom.Atom, atom2: Bio.PDB.Atom.Atom) -> bool:
    '''
    Returns 1 if two atoms share a covalent bond. Else returns 0.
    '''
    return

def covalentBondMatrix(struc: Bio.PDB.Structure) -> np.array:
    '''
    Returns the matrix that describes bonds between each pair of atoms in the structure.
    '''

    

    return

def distanceMatrix(struc:Bio.PDB.Structure) -> np.array:
    '''
    Returns the matrix that describes distances between each pair of atoms in the structure.
    '''

    return

def nonbondedDistanceWeighted(distmat: np.array()) -> np.array:
    '''
    Takes in a distance matrix and sets the array elements equal to infinity if the atoms share a covalent bond.

    Cite: Cang 2018, PLOS.
    '''