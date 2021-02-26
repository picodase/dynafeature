
'''
IMPORTS
'''

import Bio.PDB
import numpy as np

import structureutils

'''
CONSTANTS
'''



'''
FUNCTIONS
'''

class WassersteinDistance():
    def __init__(self, obj1, obj2):
        self.o1 = obj1
        self.o2 = obj2
    
    def PersistentHomology(self)-> float:
        '''
        DESCRIPTION
            Computes the Wasserstein distance between two persistent homology models of a protein structure.

        ARGUMENTS
            conf1 First conformation for comparison.
            conf2 Second filename of the second conformation for comparison.
        '''
        
        return