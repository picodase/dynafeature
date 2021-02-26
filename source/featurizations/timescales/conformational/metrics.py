
'''
IMPORTS
'''

import numpy as np

from ..conformational.confUtils import Conformation

'''
CONSTANTS
'''



'''
FUNCTIONS
'''
class Metric():
    def __init__(self, obj1, obj2):
        self.o1 = obj1
        self.o2 = obj2
        self.min = 0
        self.max = 1
        self.value = 0

class WassersteinDistance(Metric):
    def __init__(self):

class PersistentHomology(WassersteinDistance):
    '''
    DESCRIPTION
        Computes the Wasserstein distance between two persistent homology models of a protein structure.

    ARGUMENTS
        conf1 First conformation for comparison.
        conf2 Second filename of the second conformation for comparison.
    '''
    def __init__(self, conf1:Conformation, conf2:Conformation):
        assert conf1 is Conformation
        assert conf2 is Conformation
        self.o1 = conf1
        self.o2 = conf2


